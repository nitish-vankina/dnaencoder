//This is my project for the hack.club stardance engineering challenge. It implements a simple DNA encoding scheme and a custom file format for storing binary data as nucleotide strands. The library provides both in-memory and streaming APIs for encoding and decoding, with error handling and basic statistics reporting.
//For the last few years, I have found myself increasigly interested in the merging of biology and computing. This started when I first saw the Cortical Labs computer. This computer used human neurons as a CPU for a computer. I want to work in this field when I am older, and have started doing relevant projects. This project is one of them.
//I started this project because I believe that the current form of data storage is extremely inefficient, and that DNA storage is the future of data storage. I also wanted to learn more about how DNA storage works, and how to implement it in software. This project is a simple implementation of a DNA storage engine, which can encode and decode binary data as nucleotide strands. It also includes a command-line interface for encoding and decoding files, as well as a demo mode for testing the functionality of the library.

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
///   byte = b7 b6 | b5 b4 | b3 b2 | b1 b0
///           diad0   diad1   diad2   diad3
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
/// Returns `None` for any character outside {A, C, G, T, a, c, g, t}.
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
// Bulk encoder  (byte slice → nucleotide Vec)
// ---------------------------------------------------------------------------

/// Encode an in-memory byte slice into a nucleotide string (ASCII bytes).
///
/// Uses Rayon parallel iteration when `data.len() >= PARALLEL_THRESHOLD`.
///
/// # Example
/// ```
/// use helix_archive::encode_bytes;
/// // 'H' = 0x48 = 01 00 10 00 → CAGA
/// // 'i' = 0x69 = 01 10 10 01 → CGAC  (wait—let's just test round-trip)
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
/// [6..14]  payload_len u64  LE original byte count
/// [14..18] checksum    u32  LE XOR-folded checksum of header[0..14]
/// ```
pub fn encode_stream<R: Read, W: Write>(
    reader: R,
    writer: W,
    chunk_size: usize,
) -> HelixResult<EncodingStats> {
    let mut reader = BufReader::with_capacity(chunk_size, reader);
    let mut writer = BufWriter::with_capacity(chunk_size * 4, writer);

    // Buffer the encoded output so we can prepend the header with raw_bytes
    let mut encoded_buf: Vec<u8> = Vec::new();
    let mut raw_bytes: u64 = 0;
    let mut chunk = vec![0u8; chunk_size];

    loop {
        let n = reader.read(&mut chunk)?;
        if n == 0 { break; }
        raw_bytes += n as u64;
        encoded_buf.extend_from_slice(&encode_bytes(&chunk[..n]));
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

    // Write encoded strand
    writer.write_all(&encoded_buf)?;
    writer.flush()?;

    let nucleotide_count = encoded_buf.len() as u64;
    Ok(EncodingStats {
        raw_bytes,
        nucleotide_count,
        expansion_ratio: if raw_bytes > 0 {
            nucleotide_count as f64 / raw_bytes as f64
        } else { 0.0 },
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
    reader: R,
    writer: W,
    chunk_size: usize,
) -> HelixResult<DecodingStats> {
    // Align chunk to 4-base quads
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

    // Streaming decode loop.
    //
    // Key invariant: decode_bytes requires a multiple-of-4 length slice.
    // After read_exact(header), BufReader's internal buffer is misaligned
    // (capacity-18 bytes remain), so we carry a leftover ring across reads.
    let mut buf = vec![0u8; chunk_size];
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
            // Fast path: no copying needed
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
            // Slow path: merge leftover + new chunk, decode aligned prefix
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

    // Drain residual bases (valid strand => empty; included for robustness)
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
            got: decoded_bytes,
        });
    }

    Ok(DecodingStats { nucleotides_read, decoded_bytes })
}

// ---------------------------------------------------------------------------
// Stats structs
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct EncodingStats {
    pub raw_bytes:        u64,
    pub nucleotide_count: u64,
    pub expansion_ratio:  f64,
}

#[derive(Debug, Clone)]
pub struct DecodingStats {
    pub nucleotides_read: u64,
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

    #[test]
    fn encode_byte_known_values() {
        assert_eq!(encode_byte(0x00), *b"AAAA");
        assert_eq!(encode_byte(0xFF), *b"TTTT");
        // 0x6C = 0b01_10_11_00 → C G T A
        assert_eq!(encode_byte(0x6C), *b"CGTA");
        // 0x93 = 0b10_01_00_11 → G C A T
        assert_eq!(encode_byte(0x93), *b"GCAT");
    }

    #[test]
    fn decode_base_all_cases() {
        assert_eq!(decode_base(b'A'), Some(0b00));
        assert_eq!(decode_base(b'C'), Some(0b01));
        assert_eq!(decode_base(b'G'), Some(0b10));
        assert_eq!(decode_base(b'T'), Some(0b11));
        assert_eq!(decode_base(b'a'), Some(0b00));
        assert_eq!(decode_base(b'X'), None);
    }

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

    #[test]
    fn slice_bad_length() {
        assert!(decode_bytes(b"ACG").is_err()); // length 3 — not a multiple of 4
    }

    #[test]
    fn slice_invalid_base() {
        let err = decode_bytes(b"ACGX").unwrap_err();
        assert!(matches!(err, HelixError::InvalidBase('X', 3)));
    }

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

    #[test] fn stream_empty()       { stream_round_trip(&[]); }
    #[test] fn stream_one_byte()    { stream_round_trip(&[0xAB]); }
    #[test] fn stream_ascii()       { stream_round_trip(b"Stardance Engineering Challenge"); }
    #[test] fn stream_all_bytes()   { stream_round_trip(&(0u8..=255).collect::<Vec<_>>()); }
    #[test] fn stream_cross_chunk() { stream_round_trip(&vec![0xDE; CHUNK_SIZE + 13]); }

    #[test]
    fn stream_bad_magic() {
        let mut out = Vec::new();
        let err = decode_stream(Cursor::new(b"NOTAHELIX_____GARBAGE"), Cursor::new(&mut out), CHUNK_SIZE).unwrap_err();
        assert!(matches!(err, HelixError::BadMagic | HelixError::Io(_)));
    }
    #[test]
    fn stream_exact_chunk_size() {
        // Exactly CHUNK_SIZE bytes — exercises the BufReader alignment edge case
        // where read_exact(18) leaves (CHUNK_SIZE-18) bytes in the internal buffer
        stream_round_trip(&vec![0xAB; CHUNK_SIZE]);
    }

    #[test]
    fn stream_two_chunks() {
        stream_round_trip(&vec![0xCD; CHUNK_SIZE * 2]);

}