// =============================================================================
//  HelixArchive — Integration Tests
//  tests/integration.rs
// =============================================================================

use helix_archive::*;
use std::io::Cursor;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn encode_decode(data: &[u8]) -> Vec<u8> {
    let mut enc = Vec::new();
    encode_stream(Cursor::new(data), Cursor::new(&mut enc), CHUNK_SIZE)
        .expect("encode_stream failed");
    let mut dec = Vec::new();
    decode_stream(Cursor::new(&enc), Cursor::new(&mut dec), CHUNK_SIZE)
        .expect("decode_stream failed");
    dec
}

// ---------------------------------------------------------------------------
// Round-trip coverage
// ---------------------------------------------------------------------------

#[test]
fn roundtrip_empty() {
    assert_eq!(encode_decode(b""), b"");
}

#[test]
fn roundtrip_single_byte_every_value() {
    for b in 0u8..=255 {
        assert_eq!(encode_decode(&[b]), &[b], "failed for byte 0x{:02X}", b);
    }
}

#[test]
fn roundtrip_all_bytes_in_sequence() {
    let data: Vec<u8> = (0u8..=255).collect();
    assert_eq!(encode_decode(&data), data);
}

#[test]
fn roundtrip_text() {
    let msg = b"Stardance Engineering Challenge — HelixArchive";
    assert_eq!(encode_decode(msg), msg);
}

#[test]
fn roundtrip_binary_pattern() {
    let data: Vec<u8> = (0..1024).map(|i| ((i * 37 + 19) % 256) as u8).collect();
    assert_eq!(encode_decode(&data), data);
}

#[test]
fn roundtrip_cross_chunk_boundary() {
    let data = vec![0xAB_u8; CHUNK_SIZE + 1];
    assert_eq!(encode_decode(&data), data);
}

#[test]
fn roundtrip_exactly_chunk_size() {
    let data = vec![0xCD_u8; CHUNK_SIZE];
    assert_eq!(encode_decode(&data), data);
}

#[test]
fn roundtrip_multiple_chunks() {
    let data = vec![0xEF_u8; CHUNK_SIZE * 3 + 7];
    assert_eq!(encode_decode(&data), data);
}

// ---------------------------------------------------------------------------
// Header validation
// ---------------------------------------------------------------------------

#[test]
fn bad_magic_rejected() {
    let mut out = Vec::new();
    let err = decode_stream(
        Cursor::new(b"XXXX\x01\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00"),
        Cursor::new(&mut out),
        CHUNK_SIZE,
    ).unwrap_err();
    assert!(matches!(err, HelixError::BadMagic));
}

#[test]
fn bad_version_rejected() {
    // Build a valid-looking header with version=99
    let mut header = [0u8; 18];
    header[0..4].copy_from_slice(b"HXAR");
    header[4] = 99;
    let mut out = Vec::new();
    let err = decode_stream(Cursor::new(header.as_ref()), Cursor::new(&mut out), CHUNK_SIZE)
        .unwrap_err();
    assert!(matches!(err, HelixError::UnsupportedVersion(99)));
}

// ---------------------------------------------------------------------------
// EncodingStats correctness
// ---------------------------------------------------------------------------

#[test]
fn stats_expansion_ratio() {
    let data = vec![0x00u8; 256];
    let mut enc = Vec::new();
    let stats = encode_stream(Cursor::new(&data), Cursor::new(&mut enc), CHUNK_SIZE).unwrap();
    assert_eq!(stats.raw_bytes, 256);
    assert_eq!(stats.nucleotide_count, 1024);
    assert!((stats.expansion_ratio - 4.0).abs() < 1e-9);
}

#[test]
fn stats_gc_all_zeros() {
    // 0x00 → AAAA → GC = 0
    let mut enc = Vec::new();
    let stats = encode_stream(Cursor::new(&[0x00u8; 64]), Cursor::new(&mut enc), CHUNK_SIZE).unwrap();
    assert!((stats.gc_fraction).abs() < 1e-9);
}

#[test]
fn stats_gc_all_ff() {
    // 0xFF → TTTT → GC = 0
    let mut enc = Vec::new();
    let stats = encode_stream(Cursor::new(&[0xFFu8; 64]), Cursor::new(&mut enc), CHUNK_SIZE).unwrap();
    assert!((stats.gc_fraction).abs() < 1e-9);
}

#[test]
fn stats_gc_mixed() {
    // 0x6B = 01 10 10 11 → CGGT, GC=3/4=0.75
    let mut enc = Vec::new();
    let stats = encode_stream(Cursor::new(&[0x6Bu8; 64]), Cursor::new(&mut enc), CHUNK_SIZE).unwrap();
    assert!((stats.gc_fraction - 0.75).abs() < 1e-6);
}

// ---------------------------------------------------------------------------
// GC content
// ---------------------------------------------------------------------------

#[test]
fn gc_content_pure_at() {
    assert!((gc_content(b"ATATATATAT") - 0.0).abs() < 1e-9);
}

#[test]
fn gc_content_pure_gc() {
    assert!((gc_content(b"GCGCGCGCGC") - 1.0).abs() < 1e-9);
}

#[test]
fn gc_content_balanced() {
    assert!((gc_content(b"ACGT") - 0.5).abs() < 1e-9);
}

// ---------------------------------------------------------------------------
// Homopolymer detection
// ---------------------------------------------------------------------------

#[test]
fn homopolymer_finds_long_run() {
    let runs = find_homopolymer_runs(b"AAAAACGT", 5);
    assert_eq!(runs.len(), 1);
    assert_eq!(runs[0].base, 'A');
    assert_eq!(runs[0].length, 5);
    assert_eq!(runs[0].position, 0);
}

#[test]
fn homopolymer_respects_min_length() {
    let runs = find_homopolymer_runs(b"AAAGGG", 4);
    assert!(runs.is_empty()); // max run is 3, below threshold of 4
}

#[test]
fn homopolymer_multiple_runs() {
    let runs = find_homopolymer_runs(b"AAAACGGGG", 4);
    assert_eq!(runs.len(), 2);
}

#[test]
fn homopolymer_empty_strand() {
    assert!(find_homopolymer_runs(b"", 1).is_empty());
}

// ---------------------------------------------------------------------------
// Reverse complement
// ---------------------------------------------------------------------------

#[test]
fn reverse_complement_self_complementary() {
    assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
}

#[test]
fn reverse_complement_all_a() {
    assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
}

#[test]
fn reverse_complement_double_application() {
    let strand = b"GCTAGCTAGC";
    let rc = reverse_complement(strand);
    let rc_rc = reverse_complement(&rc);
    assert_eq!(rc_rc.as_slice(), strand.as_slice());
}

// ---------------------------------------------------------------------------
// Melting temperature
// ---------------------------------------------------------------------------

#[test]
fn melting_temp_all_gc_hotter() {
    let tm_gc = melting_temperature(b"GCGCGCGCGCGCGCGC").unwrap();
    let tm_at = melting_temperature(b"ATATATATATATATATATATAT").unwrap();
    // GC-rich strands should have higher Tm
    assert!(tm_gc > tm_at, "GC strand should be hotter: {} vs {}", tm_gc, tm_at);
}

#[test]
fn melting_temp_invalid_base() {
    assert!(melting_temperature(b"ACGNACGT").is_none());
}

#[test]
fn melting_temp_empty() {
    // Empty → Wallace rule: 0 AT, 0 GC → 0.0
    let tm = melting_temperature(b"").unwrap();
    assert!((tm - 0.0).abs() < 1e-9);
}

// ---------------------------------------------------------------------------
// FASTA export
// ---------------------------------------------------------------------------

#[test]
fn fasta_header_format() {
    let mut out = Vec::new();
    write_fasta(b"ACGTACGT", "seq1", Some("test strand"), 60, &mut out).unwrap();
    let text = std::str::from_utf8(&out).unwrap();
    assert!(text.starts_with(">seq1 test strand\n"));
    assert!(text.contains("ACGTACGT\n"));
}

#[test]
fn fasta_line_wrapping() {
    let strand = vec![b'A'; 200];
    let mut out = Vec::new();
    write_fasta(&strand, "seq", None, 60, &mut out).unwrap();
    let text = std::str::from_utf8(&out).unwrap();
    // First sequence line should be exactly 60 chars
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines[1].len(), 60);
    assert_eq!(lines[2].len(), 60);
    assert_eq!(lines[3].len(), 80);  // 200 - 60 - 60 = 80
}

#[test]
fn fasta_no_description() {
    let mut out = Vec::new();
    write_fasta(b"GCTA", "myseq", None, 60, &mut out).unwrap();
    let text = std::str::from_utf8(&out).unwrap();
    assert!(text.starts_with(">myseq\n"));
}

// ---------------------------------------------------------------------------
// analyze_strand
// ---------------------------------------------------------------------------

#[test]
fn analyze_strand_basic() {
    let strand = encode_bytes(b"Hello");
    let analysis = analyze_strand(&strand, 100);
    assert_eq!(analysis.length, strand.len());
    assert!(analysis.gc_fraction >= 0.0 && analysis.gc_fraction <= 1.0);
    // Reverse complement of reverse complement is identity
    let rc_rc = reverse_complement(&analysis.reverse_complement);
    assert_eq!(rc_rc, strand);
}
