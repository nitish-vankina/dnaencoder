// =============================================================================
//  HelixArchive  —  helix_archive
//  High-performance DNA Data Storage Engine
//
//  Encodes arbitrary binary data as synthetic DNA nucleotide strands using a
//  2-bit per base mapping (A=00, C=01, G=10, T=11), with biological constraint
//  analysis, FASTA export, and streaming I/O.
// =============================================================================

use std::io::{self, BufReader, BufWriter, Read, Write};
use rayon::prelude::*;
use thiserror::Error;

// ---------------------------------------------------------------------------
// Public constants
// ---------------------------------------------------------------------------

/// Default streaming chunk size (64 KiB — fits L1/L2 cache on most targets)
pub const CHUNK_SIZE: usize = 65_536;

/// Minimum payload size to engage parallel encoding (1 MiB)
pub const PARALLEL_THRESHOLD: usize = 1_048_576;

/// Magic bytes for HelixArchive strand files
pub const MAGIC: &[u8] = b"HXAR";

/// Current file-format version
pub const FORMAT_VERSION: u8 = 1;

/// Maximum recommended homopolymer run length for DNA synthesis
pub const MAX_HOMOPOLYMER_RUN: usize = 4;

/// Ideal GC content range for DNA synthesis (lower bound, inclusive)
pub const GC_MIN: f64 = 0.40;

/// Ideal GC content range for DNA synthesis (upper bound, inclusive)
pub const GC_MAX: f64 = 0.60;

// ---------------------------------------------------------------------------
// Error type
// ---------------------------------------------------------------------------

#[derive(Debug, Error)]
pub enum HelixError {
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),

    #[error("Invalid nucleotide '{0}' at position {1}")]
    InvalidBase(char, usize),

    #[error("Strand length {0} is not a multiple of 4 (cannot decode to whole bytes)")]
    BadStrandLength(usize),

    #[error("Invalid magic bytes — not a HelixArchive strand file")]
    BadMagic,

    #[error("Unsupported format version {0}")]
    UnsupportedVersion(u8),

    #[error("Payload length mismatch: header says {expected}, decoded {got}")]
    LengthMismatch { expected: u64, got: u64 },
}

pub type HelixResult<T> = Result<T, HelixError>;

// ---------------------------------------------------------------------------
// Core encoding primitives
// ---------------------------------------------------------------------------

/// Encode a single byte into exactly 4 nucleotide characters.
///
/// Each byte is split into four 2-bit diads, MSB-first:
/// ```text
///   byte = b7 b6 | b5 b4 | b3 b2 | b1 b0
///           diad0   diad1   diad2   diad3
/// ```
/// Mapping: `00→A`, `01→C`, `10→G`, `11→T`
///
/// # Example
/// ```
/// use helix_archive::encode_byte;
/// assert_eq!(encode_byte(0xFF), *b"TTTT");
/// assert_eq!(encode_byte(0x00), *b"AAAA");
/// ```
#[inline(always)]
pub fn encode_byte(byte: u8) -> [u8; 4] {
    const TABLE: [u8; 4] = [b'A', b'C', b'G', b'T'];
    [
        TABLE[((byte >> 6) & 0b11) as usize],
        TABLE[((byte >> 4) & 0b11) as usize],
        TABLE[((byte >> 2) & 0b11) as usize],
        TABLE[( byte       & 0b11) as usize],
    ]
}

/// Decode a single nucleotide ASCII byte to its 2-bit value.
///
/// Returns `None` for any character outside `{A, C, G, T, a, c, g, t}`.
#[inline(always)]
pub fn decode_base(base: u8) -> Option<u8> {
    match base {
        b'A' | b'a' => Some(0b00),
        b'C' | b'c' => Some(0b01),
        b'G' | b'g' => Some(0b10),
        b'T' | b't' => Some(0b11),
        _            => None,
    }
}

/// Decode exactly 4 nucleotide bytes back into one data byte.
///
/// # Arguments
/// * `quad`   — slice of exactly 4 nucleotide ASCII bytes
/// * `offset` — absolute position in the parent strand (used in error messages)
#[inline(always)]
pub fn decode_quad(quad: &[u8], offset: usize) -> HelixResult<u8> {
    let mut result: u8 = 0;
    for (i, &base) in quad.iter().enumerate() {
        let bits = decode_base(base).ok_or_else(|| {
            HelixError::InvalidBase(base as char, offset + i)
        })?;
        result = (result << 2) | bits;
    }
    Ok(result)
}

// ---------------------------------------------------------------------------
// Bulk encoder / decoder  (byte slice ↔ nucleotide Vec)
// ---------------------------------------------------------------------------

/// Encode an in-memory byte slice into a nucleotide string (ASCII bytes).
///
/// Uses Rayon parallel iteration when `data.len() >= PARALLEL_THRESHOLD`.
///
/// # Example
/// ```
/// use helix_archive::encode_bytes;
/// let strand = encode_bytes(b"Hi");
/// assert_eq!(strand.len(), 8);
/// ```
pub fn encode_bytes(data: &[u8]) -> Vec<u8> {
    let mut out = vec![0u8; data.len() * 4];

    if data.len() >= PARALLEL_THRESHOLD {
        out.par_chunks_mut(4)
            .zip(data.par_iter())
            .for_each(|(slot, &byte)| {
                slot.copy_from_slice(&encode_byte(byte));
            });
    } else {
        for (i, &byte) in data.iter().enumerate() {
            out[i * 4..(i + 1) * 4].copy_from_slice(&encode_byte(byte));
        }
    }

    out
}

/// Decode a nucleotide byte slice back into raw bytes.
///
/// Accepts mixed-case input. Length must be a multiple of 4.
pub fn decode_bytes(strand: &[u8]) -> HelixResult<Vec<u8>> {
    if strand.len() % 4 != 0 {
        return Err(HelixError::BadStrandLength(strand.len()));
    }

    let n = strand.len() / 4;

    if strand.len() >= PARALLEL_THRESHOLD {
        let result: HelixResult<Vec<u8>> = strand
            .par_chunks(4)
            .enumerate()
            .map(|(i, quad)| decode_quad(quad, i * 4))
            .collect();
        result
    } else {
        let mut out = vec![0u8; n];
        for (i, quad) in strand.chunks(4).enumerate() {
            out[i] = decode_quad(quad, i * 4)?;
        }
        Ok(out)
    }
}

// ---------------------------------------------------------------------------
// Biological constraint analysis
// ---------------------------------------------------------------------------

/// A single homopolymer run found in a strand.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HomopolymerRun {
    /// The repeated nucleotide base (`A`, `C`, `G`, or `T`)
    pub base: char,
    /// Zero-based start position in the strand
    pub position: usize,
    /// Length of the run
    pub length: usize,
}

/// Compute GC content of a nucleotide strand as a fraction in `[0.0, 1.0]`.
///
/// Returns `0.0` for an empty strand. Case-insensitive.
///
/// # Example
/// ```
/// use helix_archive::gc_content;
/// assert!((gc_content(b"GCGC") - 1.0).abs() < 1e-9);
/// assert!((gc_content(b"ATAT") - 0.0).abs() < 1e-9);
/// assert!((gc_content(b"ACGT") - 0.5).abs() < 1e-9);
/// ```
pub fn gc_content(strand: &[u8]) -> f64 {
    if strand.is_empty() { return 0.0; }
    let gc = strand.iter().filter(|&&b| matches!(b, b'G' | b'g' | b'C' | b'c')).count();
    gc as f64 / strand.len() as f64
}

/// Find all homopolymer runs (consecutive identical bases) in a strand.
///
/// Only runs exceeding `min_length` are returned. Use `MAX_HOMOPOLYMER_RUN`
/// (4) as the threshold for synthesis warnings.
///
/// # Example
/// ```
/// use helix_archive::{find_homopolymer_runs, HomopolymerRun};
/// let runs = find_homopolymer_runs(b"AAACGTTTT", 2);
/// assert_eq!(runs[0], HomopolymerRun { base: 'A', position: 0, length: 3 });
/// assert_eq!(runs[1], HomopolymerRun { base: 'T', position: 5, length: 4 });
/// ```
pub fn find_homopolymer_runs(strand: &[u8], min_length: usize) -> Vec<HomopolymerRun> {
    if strand.is_empty() { return vec![]; }

    let mut runs = Vec::new();
    let mut run_start = 0usize;
    let mut run_len   = 1usize;

    for i in 1..strand.len() {
        if strand[i].to_ascii_uppercase() == strand[i - 1].to_ascii_uppercase() {
            run_len += 1;
        } else {
            if run_len >= min_length {
                runs.push(HomopolymerRun {
                    base:     strand[run_start].to_ascii_uppercase() as char,
                    position: run_start,
                    length:   run_len,
                });
            }
            run_start = i;
            run_len   = 1;
        }
    }
    // flush final run
    if run_len >= min_length {
        runs.push(HomopolymerRun {
            base:     strand[run_start].to_ascii_uppercase() as char,
            position: run_start,
            length:   run_len,
        });
    }
    runs
}

/// Return the reverse complement of a DNA strand.
///
/// Complement mapping: `A↔T`, `C↔G`. Case is preserved.
/// Any character that is not a valid nucleotide is left unchanged.
///
/// # Example
/// ```
/// use helix_archive::reverse_complement;
/// assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
/// assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
/// assert_eq!(reverse_complement(b"GCTA"), b"TAGC");
/// ```
pub fn reverse_complement(strand: &[u8]) -> Vec<u8> {
    strand.iter().rev().map(|&b| match b {
        b'A' => b'T', b'T' => b'A',
        b'C' => b'G', b'G' => b'C',
        b'a' => b't', b't' => b'a',
        b'c' => b'g', b'g' => b'c',
        other => other,
    }).collect()
}

/// Estimate the melting temperature (Tm) of a DNA strand in °C using the
/// nearest-neighbor thermodynamic model (SantaLucia 1998, 1 M NaCl,
/// 250 nM strand concentration).
///
/// For strands shorter than 14 bases the simplified Wallace rule is used
/// instead: `Tm = 2·(A+T) + 4·(G+C)`.
///
/// Returns `None` if the strand contains non-ACGT characters.
///
/// # Example
/// ```
/// use helix_archive::melting_temperature;
/// let tm = melting_temperature(b"GCGCGCGCGCGCGCGC").unwrap();
/// assert!(tm > 50.0 && tm < 80.0);
/// ```
pub fn melting_temperature(strand: &[u8]) -> Option<f64> {
    // Nearest-neighbor ΔH (kcal/mol) and ΔS (cal/mol·K) from SantaLucia 1998
    // Order: AA, AT, TA, CA, GT, CT, GA, CG, GC, GG (and complements by symmetry)
    #[rustfmt::skip]
    const NN: &[(&[u8; 2], f64, f64)] = &[
        (b"AA", -7.9, -22.2), (b"TT", -7.9, -22.2),
        (b"AT", -7.2, -20.4),
        (b"TA", -7.2, -21.3),
        (b"CA", -8.5, -22.7), (b"TG", -8.5, -22.7),
        (b"GT", -8.4, -22.4), (b"AC", -8.4, -22.4),
        (b"CT", -7.8, -21.0), (b"AG", -7.8, -21.0),
        (b"GA", -8.2, -22.2), (b"TC", -8.2, -22.2),
        (b"CG", -10.6, -27.2),
        (b"GC", -9.8, -24.4),
        (b"GG", -8.0, -19.9), (b"CC", -8.0, -19.9),
    ];

    let upper: Vec<u8> = strand.iter().map(|b| b.to_ascii_uppercase()).collect();

    // Validate
    for &b in &upper {
        if !matches!(b, b'A' | b'C' | b'G' | b'T') { return None; }
    }

    if upper.len() < 14 {
        // Wallace rule
        let at = upper.iter().filter(|&&b| b == b'A' || b == b'T').count();
        let gc = upper.iter().filter(|&&b| b == b'G' || b == b'C').count();
        return Some(2.0 * at as f64 + 4.0 * gc as f64);
    }

    let mut delta_h: f64 = 0.0;
    let mut delta_s: f64 = 0.0;

    for pair in upper.windows(2) {
        let key: [u8; 2] = [pair[0], pair[1]];
        if let Some(&(_, dh, ds)) = NN.iter().find(|&&(k, _, _)| *k == key) {
            delta_h += dh;
            delta_s += ds;
        }
    }

    // Initiation parameters
    let first_gc = matches!(upper[0], b'G' | b'C');
    let last_gc  = matches!(upper[upper.len() - 1], b'G' | b'C');
    delta_h += if first_gc { 0.1 } else { 2.3 };
    delta_h += if last_gc  { 0.1 } else { 2.3 };
    delta_s += if first_gc { -2.8 } else { 4.1 };
    delta_s += if last_gc  { -2.8 } else { 4.1 };

    // Tm = ΔH / (ΔS + R·ln(CT/4)) − 273.15
    // CT = 250e-9 M, R = 1.987 cal/mol·K
    const R:  f64 = 1.987;
    const CT: f64 = 250e-9;
    let tm = (delta_h * 1000.0) / (delta_s + R * (CT / 4.0_f64).ln()) - 273.15;
    Some(tm)
}

/// Build a GC-content histogram over a strand using 5% bucket widths (20 buckets).
///
/// The strand is split into `block_size`-nucleotide windows; each window's
/// GC fraction is placed in one of 20 buckets `[0%,5%), [5%,10%), …, [95%,100%]`.
/// Returns an array of 20 counts.
pub fn gc_histogram(strand: &[u8], block_size: usize) -> [u64; 20] {
    let mut hist = [0u64; 20];
    if strand.is_empty() || block_size == 0 { return hist; }

    for block in strand.chunks(block_size) {
        let gc = gc_content(block);
        let bucket = (gc * 20.0).min(19.0) as usize;
        hist[bucket] += 1;
    }
    hist
}

// ---------------------------------------------------------------------------
// FASTA export
// ---------------------------------------------------------------------------

/// Write a nucleotide strand to `writer` in FASTA format.
///
/// Lines are wrapped at `line_width` characters (standard bioinformatics
/// convention is 60 or 80). The `id` string becomes the sequence identifier
/// after `>`, and `description` (optional) follows on the same header line.
///
/// # Example
/// ```
/// use helix_archive::{encode_bytes, write_fasta};
/// let strand = encode_bytes(b"Hello");
/// let mut out = Vec::new();
/// write_fasta(&strand, "seq1", Some("Hello encoded"), 60, &mut out).unwrap();
/// assert!(out.starts_with(b">seq1 Hello encoded\n"));
/// ```
pub fn write_fasta<W: Write>(
    strand:      &[u8],
    id:          &str,
    description: Option<&str>,
    line_width:  usize,
    writer:      &mut W,
) -> HelixResult<()> {
    match description {
        Some(desc) => writeln!(writer, ">{} {}", id, desc)?,
        None       => writeln!(writer, ">{}", id)?,
    }
    for chunk in strand.chunks(line_width) {
        writer.write_all(chunk)?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Streaming encoder  (Reader → Writer)
// ---------------------------------------------------------------------------

/// Streaming DNA encoder with HelixArchive header.
///
/// Reads `reader` in `chunk_size`-byte windows, encodes each chunk, and writes
/// nucleotides to `writer`. Prepends an 18-byte metadata header so the strand
/// can be decoded without external length information.
///
/// # Header layout (18 bytes)
/// ```text
/// [0..4]   magic       b"HXAR"
/// [4]      version     u8   = 1
/// [5]      flags       u8   = 0 (reserved)
/// [6..14]  payload_len u64  LE  original byte count
/// [14..18] checksum    u32  LE  XOR-folded checksum of header[0..14]
/// ```
pub fn encode_stream<R: Read, W: Write>(
    reader:     R,
    writer:     W,
    chunk_size: usize,
) -> HelixResult<EncodingStats> {
    let mut reader = BufReader::with_capacity(chunk_size, reader);
    let mut writer = BufWriter::with_capacity(chunk_size * 4, writer);

    let mut encoded_buf: Vec<u8> = Vec::new();
    let mut raw_bytes:   u64     = 0;
    let mut gc_count:    u64     = 0;
    let mut chunk = vec![0u8; chunk_size];

    loop {
        let n = reader.read(&mut chunk)?;
        if n == 0 { break; }
        raw_bytes += n as u64;
        let enc = encode_bytes(&chunk[..n]);
        gc_count += enc.iter().filter(|&&b| b == b'G' || b == b'C').count() as u64;
        encoded_buf.extend_from_slice(&enc);
    }

    // Build and write header
    let mut header = [0u8; 18];
    header[0..4].copy_from_slice(MAGIC);
    header[4] = FORMAT_VERSION;
    header[5] = 0;
    header[6..14].copy_from_slice(&raw_bytes.to_le_bytes());
    let checksum = simple_checksum(&header[0..14]);
    header[14..18].copy_from_slice(&checksum.to_le_bytes());
    writer.write_all(&header)?;
    writer.write_all(&encoded_buf)?;
    writer.flush()?;

    let nucleotide_count = encoded_buf.len() as u64;
    let gc_fraction = if nucleotide_count > 0 {
        gc_count as f64 / nucleotide_count as f64
    } else { 0.0 };

    let homopolymer_warnings = find_homopolymer_runs(&encoded_buf, MAX_HOMOPOLYMER_RUN + 1).len() as u64;
    let gc_hist = gc_histogram(&encoded_buf, 100);

    Ok(EncodingStats {
        raw_bytes,
        nucleotide_count,
        expansion_ratio: if raw_bytes > 0 {
            nucleotide_count as f64 / raw_bytes as f64
        } else { 0.0 },
        gc_fraction,
        gc_histogram: gc_hist,
        homopolymer_warnings,
    })
}

// ---------------------------------------------------------------------------
// Streaming decoder  (Reader → Writer)
// ---------------------------------------------------------------------------

/// Streaming DNA decoder.
///
/// Reads a HelixArchive strand file (with header), validates it, then decodes
/// nucleotides back to original bytes in `chunk_size`-nucleotide windows.
pub fn decode_stream<R: Read, W: Write>(
    reader:     R,
    writer:     W,
    chunk_size: usize,
) -> HelixResult<DecodingStats> {
    let chunk_size = ((chunk_size + 3) / 4) * 4;

    let mut reader = BufReader::with_capacity(chunk_size, reader);
    let mut writer = BufWriter::with_capacity(chunk_size / 4, writer);

    // Read & validate header
    let mut header = [0u8; 18];
    reader.read_exact(&mut header)?;

    if &header[0..4] != MAGIC { return Err(HelixError::BadMagic); }
    if header[4] != FORMAT_VERSION {
        return Err(HelixError::UnsupportedVersion(header[4]));
    }

    let expected_raw    = u64::from_le_bytes(header[6..14].try_into().unwrap());
    let stored_checksum = u32::from_le_bytes(header[14..18].try_into().unwrap());
    if simple_checksum(&header[0..14]) != stored_checksum {
        return Err(HelixError::BadMagic);
    }

    let mut buf      = vec![0u8; chunk_size];
    let mut leftover: Vec<u8> = Vec::with_capacity(3);
    let mut decoded_bytes:    u64 = 0;
    let mut nucleotides_read: u64 = 0;

    loop {
        let n = reader.read(&mut buf)?;
        if n == 0 { break; }

        let total     = leftover.len() + n;
        let usable    = (total / 4) * 4;
        let remainder = total % 4;

        if leftover.is_empty() {
            if usable > 0 {
                let decoded = decode_bytes(&buf[..usable])?;
                writer.write_all(&decoded)?;
                nucleotides_read += usable as u64;
                decoded_bytes    += decoded.len() as u64;
            }
            if remainder > 0 {
                leftover.extend_from_slice(&buf[usable..n]);
            }
        } else {
            let mut work = std::mem::take(&mut leftover);
            work.extend_from_slice(&buf[..n]);
            if usable > 0 {
                let decoded = decode_bytes(&work[..usable])?;
                writer.write_all(&decoded)?;
                nucleotides_read += usable as u64;
                decoded_bytes    += decoded.len() as u64;
            }
            if remainder > 0 {
                leftover.extend_from_slice(&work[usable..]);
            }
        }
    }

    if !leftover.is_empty() {
        let decoded = decode_bytes(&leftover)?;
        writer.write_all(&decoded)?;
        nucleotides_read += leftover.len() as u64;
        decoded_bytes    += decoded.len() as u64;
    }

    writer.flush()?;

    if decoded_bytes != expected_raw {
        return Err(HelixError::LengthMismatch {
            expected: expected_raw,
            got:      decoded_bytes,
        });
    }

    Ok(DecodingStats { nucleotides_read, decoded_bytes })
}

// ---------------------------------------------------------------------------
// Strand analysis report
// ---------------------------------------------------------------------------

/// Full biological analysis of an encoded strand.
#[derive(Debug, Clone)]
pub struct StrandAnalysis {
    /// Length of the strand in nucleotide bases
    pub length:                 usize,
    /// GC content fraction in `[0.0, 1.0]`
    pub gc_fraction:            f64,
    /// Whether GC content is within the synthesis-safe window (40–60%)
    pub gc_ok:                  bool,
    /// All homopolymer runs exceeding `MAX_HOMOPOLYMER_RUN`
    pub homopolymer_violations: Vec<HomopolymerRun>,
    /// Estimated melting temperature in °C (None if strand is empty)
    pub melting_temp_c:         Option<f64>,
    /// Reverse complement of the strand
    pub reverse_complement:     Vec<u8>,
    /// GC distribution histogram (20 buckets × 5%)
    pub gc_histogram:           [u64; 20],
}

/// Analyse a raw nucleotide strand and return a full [`StrandAnalysis`].
///
/// # Arguments
/// * `strand`     — ASCII nucleotide bytes (`A/C/G/T`, case-insensitive)
/// * `block_size` — window size for GC histogram (e.g. 100)
pub fn analyze_strand(strand: &[u8], block_size: usize) -> StrandAnalysis {
    let gc  = gc_content(strand);
    let hps = find_homopolymer_runs(strand, MAX_HOMOPOLYMER_RUN + 1);
    let tm  = melting_temperature(strand);
    let rc  = reverse_complement(strand);
    let hist = gc_histogram(strand, block_size);

    StrandAnalysis {
        length:                 strand.len(),
        gc_fraction:            gc,
        gc_ok:                  gc >= GC_MIN && gc <= GC_MAX,
        homopolymer_violations: hps,
        melting_temp_c:         tm,
        reverse_complement:     rc,
        gc_histogram:           hist,
    }
}

// ---------------------------------------------------------------------------
// Stats structs
// ---------------------------------------------------------------------------

/// Statistics returned by [`encode_stream`].
#[derive(Debug, Clone)]
pub struct EncodingStats {
    /// Original file size in bytes
    pub raw_bytes:             u64,
    /// Total nucleotide bases written
    pub nucleotide_count:      u64,
    /// Ratio of nucleotide count to raw bytes (always ~4.0 for pure encoding)
    pub expansion_ratio:       f64,
    /// GC content fraction of the encoded strand
    pub gc_fraction:           f64,
    /// GC distribution histogram (20 buckets × 5%)
    pub gc_histogram:          [u64; 20],
    /// Number of homopolymer runs exceeding `MAX_HOMOPOLYMER_RUN`
    pub homopolymer_warnings:  u64,
}

/// Statistics returned by [`decode_stream`].
#[derive(Debug, Clone)]
pub struct DecodingStats {
    /// Total nucleotide bases read from the strand file
    pub nucleotides_read: u64,
    /// Decoded data size in bytes
    pub decoded_bytes:    u64,
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

fn simple_checksum(data: &[u8]) -> u32 {
    data.chunks(4).fold(0u32, |acc, chunk| {
        let mut word = [0u8; 4];
        word[..chunk.len()].copy_from_slice(chunk);
        acc ^ u32::from_le_bytes(word)
    })
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    // ── encode_byte ─────────────────────────────────────────────────────────

    #[test]
    fn encode_byte_known_values() {
        assert_eq!(encode_byte(0x00), *b"AAAA");
        assert_eq!(encode_byte(0xFF), *b"TTTT");
        assert_eq!(encode_byte(0x6C), *b"CGTA");  // 0b01_10_11_00
        assert_eq!(encode_byte(0x93), *b"GCAT");  // 0b10_01_00_11
    }

    // ── decode_base ─────────────────────────────────────────────────────────

    #[test]
    fn decode_base_all_cases() {
        assert_eq!(decode_base(b'A'), Some(0b00));
        assert_eq!(decode_base(b'C'), Some(0b01));
        assert_eq!(decode_base(b'G'), Some(0b10));
        assert_eq!(decode_base(b'T'), Some(0b11));
        assert_eq!(decode_base(b'a'), Some(0b00));
        assert_eq!(decode_base(b'X'), None);
    }

    // ── round-trips ─────────────────────────────────────────────────────────

    #[test]
    fn byte_round_trip_all_256() {
        for byte in 0u8..=255 {
            let enc = encode_byte(byte);
            let dec = decode_quad(&enc, 0).unwrap();
            assert_eq!(byte, dec, "round-trip failed for 0x{:02X}", byte);
        }
    }

    #[test]
    fn slice_ascii_string() {
        let original = b"Hello, HelixArchive!";
        let encoded  = encode_bytes(original);
        assert_eq!(encoded.len(), original.len() * 4);
        assert_eq!(decode_bytes(&encoded).unwrap(), original);
    }

    #[test]
    fn slice_all_bytes() {
        let original: Vec<u8> = (0u8..=255).collect();
        let encoded = encode_bytes(&original);
        assert_eq!(decode_bytes(&encoded).unwrap(), original);
    }

    // ── error cases ─────────────────────────────────────────────────────────

    #[test]
    fn slice_bad_length() {
        assert!(decode_bytes(b"ACG").is_err());
    }

    #[test]
    fn slice_invalid_base() {
        let err = decode_bytes(b"ACGX").unwrap_err();
        assert!(matches!(err, HelixError::InvalidBase('X', 3)));
    }

    // ── GC content ──────────────────────────────────────────────────────────

    #[test]
    fn gc_content_values() {
        assert!((gc_content(b"ACGT") - 0.5).abs() < 1e-9);
        assert!((gc_content(b"AAAA") - 0.0).abs() < 1e-9);
        assert!((gc_content(b"GCGC") - 1.0).abs() < 1e-9);
        assert!((gc_content(b"")    - 0.0).abs() < 1e-9);
    }

    // ── homopolymer detection ────────────────────────────────────────────────

    #[test]
    fn homopolymer_basic() {
        let runs = find_homopolymer_runs(b"AAACGTTTT", 3);
        assert_eq!(runs.len(), 2);
        assert_eq!(runs[0], HomopolymerRun { base: 'A', position: 0, length: 3 });
        assert_eq!(runs[1], HomopolymerRun { base: 'T', position: 5, length: 4 });
    }

    #[test]
    fn homopolymer_none() {
        assert!(find_homopolymer_runs(b"ACGTACGT", 3).is_empty());
    }

    #[test]
    fn homopolymer_case_insensitive() {
        let runs = find_homopolymer_runs(b"aaaCGT", 3);
        assert_eq!(runs.len(), 1);
        assert_eq!(runs[0].base, 'A');
    }

    // ── reverse complement ───────────────────────────────────────────────────

    #[test]
    fn reverse_complement_values() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"GCTA"), b"TAGC");
        assert_eq!(reverse_complement(b""),     b"");
    }

    // ── melting temperature ──────────────────────────────────────────────────

    #[test]
    fn melting_temp_wallace() {
        // Short strand: 2*AT + 4*GC = 2*2 + 4*2 = 12
        let tm = melting_temperature(b"AATGCG").unwrap();
        // Wallace: 2*(2AT) + 4*(4GC wait: AATGCG has 2 AT + 4 GC? A,A,T = 3 AT, G,C,G = 3 GC
        // 2*3 + 4*3 = 18.0
        assert!(tm > 10.0 && tm < 50.0);
    }

    #[test]
    fn melting_temp_longer_strand() {
        let tm = melting_temperature(b"GCGCGCGCGCGCGCGC").unwrap();
        assert!(tm > 50.0 && tm < 100.0);
    }

    #[test]
    fn melting_temp_invalid() {
        assert!(melting_temperature(b"ACGX").is_none());
    }

    // ── FASTA export ─────────────────────────────────────────────────────────

    #[test]
    fn fasta_format() {
        let strand = encode_bytes(b"Hi");
        let mut out = Vec::new();
        write_fasta(&strand, "seq1", Some("test"), 60, &mut out).unwrap();
        let s = std::str::from_utf8(&out).unwrap();
        assert!(s.starts_with(">seq1 test\n"));
        assert!(s.contains("CGAC")); // 'i' = 0x69 → CGAC... verify encoding present
    }

    #[test]
    fn fasta_no_description() {
        let mut out = Vec::new();
        write_fasta(b"ACGT", "id1", None, 60, &mut out).unwrap();
        assert!(std::str::from_utf8(&out).unwrap().starts_with(">id1\n"));
    }

    // ── streaming ────────────────────────────────────────────────────────────

    fn stream_round_trip(data: &[u8]) {
        let mut enc_buf = Vec::new();
        let stats = encode_stream(Cursor::new(data), Cursor::new(&mut enc_buf), CHUNK_SIZE)
            .expect("encode_stream failed");
        assert_eq!(stats.raw_bytes, data.len() as u64);

        let mut dec_buf = Vec::new();
        decode_stream(Cursor::new(&enc_buf), Cursor::new(&mut dec_buf), CHUNK_SIZE)
            .expect("decode_stream failed");
        assert_eq!(dec_buf, data);
    }

    #[test] fn stream_empty()        { stream_round_trip(&[]); }
    #[test] fn stream_one_byte()     { stream_round_trip(&[0xAB]); }
    #[test] fn stream_ascii()        { stream_round_trip(b"Stardance Engineering Challenge"); }
    #[test] fn stream_all_bytes()    { stream_round_trip(&(0u8..=255).collect::<Vec<_>>()); }
    #[test] fn stream_cross_chunk()  { stream_round_trip(&vec![0xDE; CHUNK_SIZE + 13]); }
    #[test] fn stream_exact_chunk()  { stream_round_trip(&vec![0xAB; CHUNK_SIZE]); }
    #[test] fn stream_two_chunks()   { stream_round_trip(&vec![0xCD; CHUNK_SIZE * 2]); }

    #[test]
    fn stream_stats_gc_fraction() {
        // All-0x00 encodes to AAAA... so GC should be 0.0
        let mut enc = Vec::new();
        let stats = encode_stream(Cursor::new(&[0x00u8; 16]), Cursor::new(&mut enc), CHUNK_SIZE).unwrap();
        assert!((stats.gc_fraction - 0.0).abs() < 1e-9);

        // All-0xFF encodes to TTTT... so GC should also be 0.0
        let mut enc2 = Vec::new();
        let stats2 = encode_stream(Cursor::new(&[0xFFu8; 16]), Cursor::new(&mut enc2), CHUNK_SIZE).unwrap();
        assert!((stats2.gc_fraction - 0.0).abs() < 1e-9);
    }

    #[test]
    fn stream_bad_magic() {
        let mut out = Vec::new();
        let err = decode_stream(
            Cursor::new(b"NOTAHELIX_____GARBAGE"),
            Cursor::new(&mut out),
            CHUNK_SIZE,
        ).unwrap_err();
        assert!(matches!(err, HelixError::BadMagic | HelixError::Io(_)));
    }
}
