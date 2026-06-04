# HelixArchive

**High-performance DNA Data Storage Engine** — Stardance Engineering Challenge

HelixArchive encodes arbitrary binary data as synthetic DNA nucleotide strands,
and decodes them back without loss. It combines a fast streaming I/O engine with
biological constraint analysis, making it useful both as a learning tool and as a
foundation for real DNA data storage research.

---

## Why DNA Storage?

Conventional storage media (HDDs, flash) degrade in decades. DNA can survive
thousands of years, stores ~215 petabytes per gram, and requires no power to
maintain. Organisations including Microsoft and the Broad Institute are actively
researching DNA as a long-term archival medium.

HelixArchive was built to explore that intersection of biology and computing —
inspired by work like [Cortical Labs](https://corticalabs.com/) and the broader
field of biological computing.

---

## Encoding Scheme

Each byte is split into four 2-bit *diads*, encoded MSB-first:

```
byte = b7 b6 | b5 b4 | b3 b2 | b1 b0
        diad0   diad1   diad2   diad3

00 → A (Adenine)
01 → C (Cytosine)
10 → G (Guanine)
11 → T (Thymine)
```

**Example:** `'R'` (0x52 = `01 01 00 10`) → `CCAG`

Every byte becomes exactly 4 bases, giving a 4× expansion ratio. The encoded
strand is wrapped in an 18-byte binary header for streaming decode without
external metadata.

### File Format Header (18 bytes)

| Offset | Size | Field        | Description                          |
|--------|------|--------------|--------------------------------------|
| 0      | 4    | magic        | `HXAR`                               |
| 4      | 1    | version      | `1`                                  |
| 5      | 1    | flags        | Reserved, `0`                        |
| 6      | 8    | payload_len  | Original byte count (u64 LE)         |
| 14     | 4    | checksum     | XOR-folded checksum of bytes 0–13    |

---

## Biological Constraints

Real DNA synthesis has physical limits. HelixArchive analyses every strand for:

| Constraint | Threshold | Why it matters |
|---|---|---|
| **GC content** | 40–60% | Low/high GC causes mis-folding and synthesis failure |
| **Homopolymer runs** | ≤ 4 identical bases | Longer runs cause polymerase slippage errors |
| **Melting temperature** | Reported in °C | Affects hybridisation and sequencing accuracy |

---

## Installation

**Prerequisites:** Rust 1.70+ (`rustup.rs`)

```bash
git clone https://github.com/yourname/helix_archive
cd helix_archive
cargo build --release
```

The binary will be at `target/release/helix_cli`.

---

## Usage

### Encode a file

```bash
helix_cli encode myfile.pdf myfile.strand
```

Output:
```
  [ENCODE]  myfile.pdf → myfile.strand
  ──────────────────────────────────────────────────
  ✓ Raw size          : 1.23 MiB
  ✓ Strand length     : 5123456 bases
  ✓ Expansion ratio   : 4.0×
  ✓ Throughput        : 312.4 MiB/s
  ✓ Wall time         : 0.004 s
  ──────────────────────────────────────────────────
  Biological analysis:
  ✓ GC content        : [████████████░░░░░░░░] 52.3%
    ✓ Within synthesis range (40–60%)
  ✓ Homopolymer runs  : none exceeding 4 bases
```

### Encode + export FASTA

```bash
helix_cli encode myfile.pdf myfile.strand --fasta myfile.fasta
```

FASTA output is compatible with standard bioinformatics tools (BLAST, Benchling,
Geneious, etc.).

### Decode a strand file

```bash
helix_cli decode myfile.strand recovered.pdf
```

### Analyse a strand file

```bash
helix_cli analyze myfile.strand
```

Output includes GC content, homopolymer violations, estimated melting temperature,
and a GC distribution histogram.

### Interactive demo

```bash
helix_cli demo
```

Walks through the encoding scheme, runs a streaming benchmark, and demonstrates
error handling.

---

## Library API

```rust
use helix_archive::*;

// Single byte
let bases = encode_byte(0x52);          // → [b'C', b'C', b'A', b'G']
let byte  = decode_quad(&bases, 0)?;    // → 0x52

// Byte slice
let strand    = encode_bytes(b"Hello");
let recovered = decode_bytes(&strand)?;

// Streaming (large files)
let stats = encode_stream(reader, writer, CHUNK_SIZE)?;
println!("GC: {:.1}%", stats.gc_fraction * 100.0);

// Biological analysis
let gc   = gc_content(&strand);
let runs = find_homopolymer_runs(&strand, 5);  // runs longer than 5
let tm   = melting_temperature(&strand);       // nearest-neighbor model
let rc   = reverse_complement(&strand);

// Full analysis report
let analysis = analyze_strand(&strand, 100);

// FASTA export
write_fasta(&strand, "seq1", Some("my sequence"), 60, &mut writer)?;
```

---

## Architecture

```
helix_archive/
├── src/
│   ├── lib.rs          Core library — encoding, decoding, bio-analysis
│   └── main.rs         CLI — encode / decode / analyze / demo
├── tests/
│   └── integration.rs  End-to-end integration tests
├── benches/
│   └── throughput.rs   Criterion benchmarks
└── Cargo.toml
```

**Key design decisions:**

- **Rayon parallelism** kicks in above 1 MiB for both encode and decode, saturating
  all available cores automatically.
- **Streaming I/O** with `BufReader`/`BufWriter` keeps memory usage constant
  regardless of file size.
- **Leftover ring buffer** in the decoder handles BufReader alignment edge cases
  after the 18-byte header read.
- **Nearest-neighbor Tm model** (SantaLucia 1998) is the gold standard used in
  actual oligonucleotide synthesis pipelines.

---

## Running Tests

```bash
cargo test                 # unit + integration tests
cargo test --release       # with optimisations
cargo bench                # Criterion throughput benchmarks
```

---

## Biological Realism & Future Work

This implementation covers the software layer of DNA storage. A production system
would also need:

- **Error-correcting codes** (Reed-Solomon or fountain codes) for sequencing noise
- **GC-balancing encoding** (e.g. Goldman et al. base-3 scheme) to guarantee
  synthesis-safe strands regardless of input
- **Primer flanking** to allow PCR amplification of stored strands
- **Index addressing** so individual files can be retrieved from a pool

These are the next planned additions to HelixArchive.

---

## References

- Church, G. et al. (2012). *Next-generation digital information storage in DNA*. Science.
- Goldman, N. et al. (2013). *Towards practical, high-capacity, low-maintenance information storage in synthesized DNA*. Nature.
- SantaLucia, J. (1998). *A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics*. PNAS.
- Organick, L. et al. (2018). *Random access in large-scale DNA data storage*. Nature Biotechnology.

---

## License

MIT
