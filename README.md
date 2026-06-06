# 🧬 dnatools

> Encode any text into a real DNA sequence — and decode it back.

dnatools is a browser-based tool that converts text into DNA using 2-bit encoding. No backend, no dependencies, just two files.

[→ Live Demo](#) *(add your GitHub Pages link)*

---

## What It Does

DNA can store data. Each base (A, C, G, T) holds 2 bits, so any text can be mapped to a real nucleotide sequence and back. This project makes that process visual and interactive.

**Four tools in one:**

| Tab | What it does |
|---|---|
| **Encode** | Converts text → DNA sequence with live GC% stats |
| **Decode** | Paste a sequence and get the original text back |
| **Analyze** | GC content, melting temperature, homopolymer runs, reverse complement |
| **Trace** | Color-coded view of how each character becomes four bases |

---

## How the Encoding Works

Every character becomes a byte. Every byte splits into four 2-bit chunks. Each chunk maps to a base.

```
"H"  →  72  →  01 00 10 00  →  C A G A
"i"  →  105 →  01 10 10 01  →  C G G C
```

| Bits | Base |
|------|------|
| `00` | A (Adenine) |
| `01` | C (Cytosine) |
| `10` | G (Guanine) |
| `11` | T (Thymine) |

One character always produces exactly four bases. The encoding is fully reversible with no extra metadata.

---

## Biological Metrics

The analyzer computes real synthesis-relevant metrics:

- **GC content** — 40–60% is ideal for synthesis reliability
- **Melting temperature (Tm)** — SantaLucia 1998 nearest-neighbor model for strands ≥14 bases, Wallace rule for shorter ones
- **Homopolymer runs** — long repeats of the same base that cause synthesis errors
- **Reverse complement** — the complementary strand (A↔T, C↔G), reversed

These are the same checks that real DNA synthesis companies like Twist Bioscience and IDT run before manufacturing a strand.

---

## Why DNA Storage?

DNA is the densest storage medium ever discovered. One gram can theoretically hold ~215 petabytes. Microsoft, the Wyss Institute at Harvard, and others are actively building DNA storage systems. This project implements the same base encoding layer those systems use, with a working biological analysis engine on top.

---

## Getting Started

No install needed.

```bash
git clone https://github.com/yourusername/dnatools
cd dnatools
open index.html
```

---

## File Structure

```
dnatools/
├── index.html   # UI — tabs, inputs, visualizations
└── engine.js    # All logic — works in browser or Node.js
```

---

## Engine API

```js
encode("Hello")           // → DNA strand
decode("CAGACGAC...")     // → original text

gcContent(strand)         // → 0.0–1.0
analyzeStrand(strand)     // → { gcFraction, tm, hpRuns, revComp, baseCounts }
meltingTemp(strand)       // → °C
reverseComplement(strand) // → complementary strand
```

---

## Built With

- Vanilla HTML, CSS, JavaScript
- [Space Grotesk](https://fonts.google.com/specimen/Space+Grotesk) (Google Fonts)
- Zero runtime dependencies

---

MIT License