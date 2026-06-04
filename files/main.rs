// =============================================================================
//  HelixArchive CLI  — helix_cli
//  Stardance Engineering Challenge
//
//  Usage:
//    helix_cli encode <input> <output.strand> [--fasta <output.fasta>]
//    helix_cli decode <input.strand> <output>
//    helix_cli analyze <input.strand>
//    helix_cli demo
// =============================================================================

use helix_archive::*;
use std::{env, fs::File, io::Cursor, time::Instant};

// ---------------------------------------------------------------------------
// Formatting helpers
// ---------------------------------------------------------------------------

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

fn pct(f: f64) -> String { format!("{:.1}%", f * 100.0) }

fn gc_bar(fraction: f64, width: usize) -> String {
    let filled = (fraction * width as f64).round() as usize;
    let empty  = width.saturating_sub(filled);
    format!("[{}{}] {}", "█".repeat(filled), "░".repeat(empty), pct(fraction))
}

fn print_separator(ch: char, width: usize) {
    println!("  {}", ch.to_string().repeat(width));
}

// ---------------------------------------------------------------------------
// Commands
// ---------------------------------------------------------------------------

fn cmd_encode(input: &str, output: &str, fasta_out: Option<&str>) -> HelixResult<()> {
    println!("\n  [ENCODE]  {} → {}", input, output);
    let src = File::open(input)?;
    let dst = File::create(output)?;
    let t0  = Instant::now();
    let stats = encode_stream(src, dst, CHUNK_SIZE)?;
    let elapsed = t0.elapsed().as_secs_f64();

    print_separator('─', 50);
    println!("  ✓ Raw size          : {}", human_bytes(stats.raw_bytes));
    println!("  ✓ Strand length     : {} bases", stats.nucleotide_count);
    println!("  ✓ Expansion ratio   : {:.1}×", stats.expansion_ratio);
    println!("  ✓ Throughput        : {:.1} MiB/s", throughput_mbs(stats.raw_bytes, elapsed));
    println!("  ✓ Wall time         : {:.3} s", elapsed);
    print_separator('─', 50);
    println!("  Biological analysis:");
    println!("  ✓ GC content        : {}", gc_bar(stats.gc_fraction, 20));
    let gc_status = if stats.gc_fraction >= GC_MIN && stats.gc_fraction <= GC_MAX {
        "✓ Within synthesis range (40–60%)"
    } else {
        "⚠ Outside ideal synthesis range (40–60%)"
    };
    println!("    {}", gc_status);
    if stats.homopolymer_warnings > 0 {
        println!("  ⚠ Homopolymer runs  : {} run(s) > {} bases (synthesis risk)",
            stats.homopolymer_warnings, MAX_HOMOPOLYMER_RUN);
    } else {
        println!("  ✓ Homopolymer runs  : none exceeding {} bases", MAX_HOMOPOLYMER_RUN);
    }
    print_separator('─', 50);

    // Optional FASTA export
    if let Some(fasta_path) = fasta_out {
        let strand_file = std::fs::read(output)?;
        // Strip the 18-byte header to get raw nucleotides
        let nucleotides = &strand_file[18..];
        let mut f = File::create(fasta_path)?;
        let desc = format!("HelixArchive encoded {} ({} bases)", input, stats.nucleotide_count);
        write_fasta(nucleotides, "helix_seq_1", Some(&desc), 60, &mut f)?;
        println!("  ✓ FASTA exported    : {}", fasta_path);
    }

    Ok(())
}

fn cmd_decode(input: &str, output: &str) -> HelixResult<()> {
    println!("\n  [DECODE]  {} → {}", input, output);
    let src = File::open(input)?;
    let dst = File::create(output)?;
    let t0  = Instant::now();
    let stats = decode_stream(src, dst, CHUNK_SIZE)?;
    let elapsed = t0.elapsed().as_secs_f64();

    print_separator('─', 50);
    println!("  ✓ Bases read        : {}", stats.nucleotides_read);
    println!("  ✓ Decoded size      : {}", human_bytes(stats.decoded_bytes));
    println!("  ✓ Throughput        : {:.1} MiB/s", throughput_mbs(stats.decoded_bytes, elapsed));
    println!("  ✓ Wall time         : {:.3} s", elapsed);
    print_separator('─', 50);
    Ok(())
}

fn cmd_analyze(input: &str) -> HelixResult<()> {
    println!("\n  [ANALYZE] {}", input);

    // Read the strand file, strip header
    let raw = std::fs::read(input)?;
    if raw.len() < 18 {
        return Err(HelixError::BadMagic);
    }
    // Validate magic
    if &raw[0..4] != MAGIC { return Err(HelixError::BadMagic); }

    let expected_raw = u64::from_le_bytes(raw[6..14].try_into().unwrap());
    let nucleotides  = &raw[18..];

    let analysis = analyze_strand(nucleotides, 100);

    print_separator('─', 50);
    println!("  Strand length       : {} bases ({} bytes encoded)",
        analysis.length, human_bytes(expected_raw));
    println!("  GC content          : {}", gc_bar(analysis.gc_fraction, 20));
    println!("  Synthesis window    : {}",
        if analysis.gc_ok { "✓ 40–60% (safe)" } else { "⚠ Outside 40–60% (risk)" });

    // Melting temperature
    // For very long strands use a representative 100-base window
    let sample = &nucleotides[..nucleotides.len().min(100)];
    let tm_label = if let Some(tm) = melting_temperature(sample) {
        format!("{:.1} °C (first 100 bases)", tm)
    } else {
        "N/A".to_string()
    };
    println!("  Melting temp (Tm)   : {}", tm_label);

    // Homopolymer runs
    let violations = &analysis.homopolymer_violations;
    if violations.is_empty() {
        println!("  Homopolymer runs    : ✓ none > {} bases", MAX_HOMOPOLYMER_RUN);
    } else {
        println!("  Homopolymer runs    : ⚠ {} violation(s) > {} bases",
            violations.len(), MAX_HOMOPOLYMER_RUN);
        for (i, run) in violations.iter().take(5).enumerate() {
            println!("    [{}] '{}' ×{} at position {}", i + 1, run.base, run.length, run.position);
        }
        if violations.len() > 5 {
            println!("    … and {} more", violations.len() - 5);
        }
    }

    // GC histogram
    print_separator('─', 50);
    println!("  GC% distribution (per 100-base window):\n");
    let max_count = *analysis.gc_histogram.iter().max().unwrap_or(&1).max(&1);
    for (i, &count) in analysis.gc_histogram.iter().enumerate() {
        let label = format!("{:3}–{:3}%", i * 5, (i + 1) * 5);
        let bar_len = (count as f64 / max_count as f64 * 30.0).round() as usize;
        println!("  {} │{}{} {}",
            label,
            "█".repeat(bar_len),
            " ".repeat(30 - bar_len),
            count);
    }
    print_separator('─', 50);

    Ok(())
}

fn cmd_demo() -> HelixResult<()> {
    println!("\n┌─────────────────────────────────────────────────────────┐");
    println!("│          HelixArchive — Interactive Demo                │");
    println!("└─────────────────────────────────────────────────────────┘");

    // ── 1. 2-bit mapping table ───────────────────────────────────────────────
    println!("\n▸ 2-bit Mapping Table");
    println!("  00 → A   (Adenine)");
    println!("  01 → C   (Cytosine)");
    println!("  10 → G   (Guanine)");
    println!("  11 → T   (Thymine)");

    // ── 2. Single byte walkthrough ───────────────────────────────────────────
    let demo_byte = b'R'; // 0x52 = 01 01 00 10
    let enc4 = encode_byte(demo_byte);
    println!("\n▸ Single-byte example: '{}' (0x{:02X} = {:08b})", demo_byte as char, demo_byte, demo_byte);
    println!("  Diads:  {:02b} | {:02b} | {:02b} | {:02b}",
        (demo_byte >> 6) & 3, (demo_byte >> 4) & 3,
        (demo_byte >> 2) & 3,  demo_byte        & 3);
    println!("  Strand: {}", std::str::from_utf8(&enc4).unwrap());
    let back = decode_quad(&enc4, 0)?;
    println!("  Decoded back: '{}' (0x{:02X})  ✓", back as char, back);

    // ── 3. String encode/decode ──────────────────────────────────────────────
    let msg = b"HELIX";
    println!("\n▸ String encode: \"{}\"", std::str::from_utf8(msg).unwrap());
    let strand = encode_bytes(msg);
    let strand_str = std::str::from_utf8(&strand).unwrap();
    println!("  DNA strand : {}", strand_str);
    let rc = reverse_complement(&strand);
    println!("  Rev. comp  : {}", std::str::from_utf8(&rc).unwrap());
    let recovered = decode_bytes(&strand)?;
    println!("  Decoded    : \"{}\"  ✓", std::str::from_utf8(&recovered).unwrap());

    // ── 4. Biological analysis ───────────────────────────────────────────────
    println!("\n▸ Biological analysis of \"HELIX\" strand");
    let gc = gc_content(&strand);
    println!("  GC content  : {}", gc_bar(gc, 20));
    let runs = find_homopolymer_runs(&strand, MAX_HOMOPOLYMER_RUN + 1);
    if runs.is_empty() {
        println!("  Homopolymer : ✓ no runs > {} bases", MAX_HOMOPOLYMER_RUN);
    } else {
        for r in &runs {
            println!("  Homopolymer : ⚠ '{}' ×{} at pos {}", r.base, r.length, r.position);
        }
    }
    if let Some(tm) = melting_temperature(&strand) {
        println!("  Melting Tm  : {:.1} °C", tm);
    }

    // ── 5. Streaming benchmark ───────────────────────────────────────────────
    let sizes: &[usize] = &[1_024, 65_536, 1_048_576];
    println!("\n▸ Streaming benchmark");
    println!("  {:>10}  {:>14}  {:>12}  {:>10}  {:>8}",
        "Input", "Bases Out", "Throughput", "GC%", "HP warns");
    println!("  {}", "─".repeat(62));

    for &sz in sizes {
        let payload: Vec<u8> = (0..sz).map(|i| (i * 7 + 13) as u8).collect();
        let mut enc_buf = Vec::with_capacity(sz * 4 + 32);

        let t0 = Instant::now();
        let stats = encode_stream(Cursor::new(&payload), Cursor::new(&mut enc_buf), CHUNK_SIZE)?;
        let elapsed = t0.elapsed().as_secs_f64().max(1e-9);

        let mut dec_buf = Vec::with_capacity(sz);
        decode_stream(Cursor::new(&enc_buf), Cursor::new(&mut dec_buf), CHUNK_SIZE)?;
        assert_eq!(dec_buf, payload, "round-trip mismatch at size {}", sz);

        println!("  {:>10}  {:>14}  {:>9.1} MiB/s  {:>9}  {:>8}",
            human_bytes(stats.raw_bytes),
            format!("{} bases", stats.nucleotide_count),
            throughput_mbs(stats.raw_bytes, elapsed),
            pct(stats.gc_fraction),
            stats.homopolymer_warnings);
    }
    println!("\n  All round-trips verified ✓");

    // ── 6. Error handling ────────────────────────────────────────────────────
    println!("\n▸ Error handling");
    let bad_strand = b"ACGXACGT";
    match decode_bytes(bad_strand) {
        Err(HelixError::InvalidBase(ch, pos)) =>
            println!("  InvalidBase('{ch}', {pos}) correctly caught ✓"),
        other => println!("  Unexpected: {other:?}"),
    }

    let odd_strand = b"ACGACG";
    match decode_bytes(odd_strand) {
        Err(HelixError::BadStrandLength(n)) =>
            println!("  BadStrandLength({n}) correctly caught ✓"),
        other => println!("  Unexpected: {other:?}"),
    }

    println!("\n  Demo complete.\n");
    Ok(())
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

fn main() {
    let args: Vec<String> = env::args().collect();

    let result = match args.get(1).map(|s| s.as_str()) {
        Some("encode") if args.len() >= 4 => {
            // Optional --fasta <path> flag
            let fasta = args.windows(2).find(|w| w[0] == "--fasta").map(|w| w[1].as_str());
            cmd_encode(&args[2], &args[3], fasta)
        }
        Some("decode")  if args.len() == 4 => cmd_decode(&args[2], &args[3]),
        Some("analyze") if args.len() == 3 => cmd_analyze(&args[2]),
        Some("demo")                        => cmd_demo(),
        _ => {
            eprintln!("\nHelixArchive DNA Storage Engine");
            eprintln!("Usage:");
            eprintln!("  helix_cli encode  <input_file> <output.strand> [--fasta <output.fasta>]");
            eprintln!("  helix_cli decode  <input.strand> <output_file>");
            eprintln!("  helix_cli analyze <input.strand>");
            eprintln!("  helix_cli demo\n");
            std::process::exit(1);
        }
    };

    if let Err(e) = result {
        eprintln!("\n  Error: {}", e);
        std::process::exit(1);
    }
}
