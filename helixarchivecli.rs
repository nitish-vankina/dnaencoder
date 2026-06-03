// =============================================================================
//  HelixArchive CLI  — helix_cli
//  Stardance Engineering Challenge
//
//  Usage:
//    helix_cli encode <input> <output.strand>
//    helix_cli decode <input.strand> <output>
//    helix_cli demo
// =============================================================================

use helix_archive::*;
use std::{env, fs::File, io::Cursor, time::Instant};

fn human_bytes(n: u64) -> String {
    match n {
        b if b >= 1_073_741_824 => format!("{:.2} GiB", b as f64 / 1_073_741_824.0),
        b if b >= 1_048_576     => format!("{:.2} MiB", b as f64 / 1_048_576.0),
        b if b >= 1_024         => format!("{:.2} KiB", b as f64 / 1_024.0),
        b                       => format!("{} B", b),
    }
}

fn throughput_mbs(bytes: u64, elapsed_secs: f64) -> f64 {
    (bytes as f64 / 1_048_576.0) / elapsed_secs
}

fn cmd_encode(input: &str, output: &str) -> HelixResult<()> {
    println!("  [ENCODE]  {} → {}", input, output);
    let src = File::open(input)?;
    let dst = File::create(output)?;
    let t0 = Instant::now();
    let stats = encode_stream(src, dst, CHUNK_SIZE)?;
    let elapsed = t0.elapsed().as_secs_f64();

    println!("  ✓ Raw size      : {}", human_bytes(stats.raw_bytes));
    println!("  ✓ Strand length : {} bases", stats.nucleotide_count);
    println!("  ✓ Expansion     : {:.1}×", stats.expansion_ratio);
    println!("  ✓ Throughput    : {:.1} MiB/s", throughput_mbs(stats.raw_bytes, elapsed));
    println!("  ✓ Wall time     : {:.3} s", elapsed);
    Ok(())
}

fn cmd_decode(input: &str, output: &str) -> HelixResult<()> {
    println!("  [DECODE]  {} → {}", input, output);
    let src = File::open(input)?;
    let dst = File::create(output)?;
    let t0 = Instant::now();
    let stats = decode_stream(src, dst, CHUNK_SIZE)?;
    let elapsed = t0.elapsed().as_secs_f64();

    println!("  ✓ Bases read    : {}", stats.nucleotides_read);
    println!("  ✓ Decoded size  : {}", human_bytes(stats.decoded_bytes));
    println!("  ✓ Throughput    : {:.1} MiB/s", throughput_mbs(stats.decoded_bytes, elapsed));
    println!("  ✓ Wall time     : {:.3} s", elapsed);
    Ok(())
}

fn cmd_demo() -> HelixResult<()> {
    println!("\n┌─────────────────────────────────────────────────────────┐");
    println!("│          HelixArchive — Interactive Demo                │");
    println!("└─────────────────────────────────────────────────────────┘\n");

    // ── 1. Single byte walkthrough ──────────────────────────────────────────
    println!("▸ 2-bit Mapping Table");
    println!("  00 → A   (Adenine)");
    println!("  01 → C   (Cytosine)");
    println!("  10 → G   (Guanine)");
    println!("  11 → T   (Thymine)\n");

    let demo_byte = b'R'; // 0x52 = 01 01 00 10
    let enc4 = encode_byte(demo_byte);
    println!("▸ Single-byte example: '{}' (0x{:02X} = {:08b})", demo_byte as char, demo_byte, demo_byte);
    println!("  Diads:  {:02b} | {:02b} | {:02b} | {:02b}",
        (demo_byte >> 6) & 3, (demo_byte >> 4) & 3,
        (demo_byte >> 2) & 3,  demo_byte        & 3);
    println!("  Strand: {}", std::str::from_utf8(&enc4).unwrap());
    let back = decode_quad(&enc4, 0)?;
    println!("  Decoded back: '{}' (0x{:02X})  ✓\n", back as char, back);

    // ── 2. String encode/decode ─────────────────────────────────────────────
    let msg = b"HELIX";
    println!("▸ String encode: \"{}\"", std::str::from_utf8(msg).unwrap());
    let strand = encode_bytes(msg);
    println!("  DNA strand: {}", std::str::from_utf8(&strand).unwrap());
    let recovered = decode_bytes(&strand)?;
    println!("  Decoded:    \"{}\"  ✓\n", std::str::from_utf8(&recovered).unwrap());

    // ── 3. Streaming benchmark ──────────────────────────────────────────────
    let sizes: &[usize] = &[1_024, 65_536, 1_048_576];
    println!("▸ Streaming benchmark\n  {:>10}  {:>14}  {:>12}", "Input", "Bases Out", "Throughput");
    println!("  {}", "─".repeat(42));

    for &sz in sizes {
        let payload: Vec<u8> = (0..sz).map(|i| (i * 7 + 13) as u8).collect();
        let mut enc_buf = Vec::with_capacity(sz * 4 + 32);

        let t0 = Instant::now();
        let stats = encode_stream(Cursor::new(&payload), Cursor::new(&mut enc_buf), CHUNK_SIZE)?;
        let elapsed = t0.elapsed().as_secs_f64().max(1e-9);

        let mut dec_buf = Vec::with_capacity(sz);
        decode_stream(Cursor::new(&enc_buf), Cursor::new(&mut dec_buf), CHUNK_SIZE)?;
        assert_eq!(dec_buf, payload, "round-trip mismatch at size {}", sz);

        println!("  {:>10}  {:>14}  {:>9.1} MiB/s",
            human_bytes(stats.raw_bytes),
            format!("{} bases", stats.nucleotide_count),
            throughput_mbs(stats.raw_bytes, elapsed));
    }

    println!("\n  All round-trips verified ✓");

    // ── 4. Error handling ───────────────────────────────────────────────────
    println!("\n▸ Error handling");
    let bad_strand = b"ACGXACGT";
    match decode_bytes(bad_strand) {
        Err(HelixError::InvalidBase(ch, pos)) =>
            println!("  InvalidBase('{ch}', {pos}) correctly caught ✓"),
        other => println!("  Unexpected: {other:?}"),
    }

    let odd_strand = b"ACGACG"; // length 6 — not multiple of 4
    match decode_bytes(odd_strand) {
        Err(HelixError::BadStrandLength(n)) =>
            println!("  BadStrandLength({n}) correctly caught ✓"),
        other => println!("  Unexpected: {other:?}"),
    }

    println!("\n  Demo complete.\n");
    Ok(())
}

fn main() {
    let args: Vec<String> = env::args().collect();

    let result = match args.get(1).map(|s| s.as_str()) {
        Some("encode") if args.len() == 4 => cmd_encode(&args[2], &args[3]),
        Some("decode") if args.len() == 4 => cmd_decode(&args[2], &args[3]),
        Some("demo")                       => cmd_demo(),
        _ => {
            eprintln!("\nHelixArchive DNA Storage Engine");
            eprintln!("Usage:");
            eprintln!("  helix_cli encode <input_file> <output.strand>");
            eprintln!("  helix_cli decode <input.strand> <output_file>");
            eprintln!("  helix_cli demo\n");
            std::process::exit(1);
        }
    };

    if let Err(e) = result {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}