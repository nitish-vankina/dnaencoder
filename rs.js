/**
 * rs.js — Reed-Solomon codec over GF(256)
 *
 * Primitive polynomial: x^8 + x^4 + x^3 + x^2 + 1  (0x11d)
 * Generator:            α = 0x02 (root of the primitive poly)
 *
 * Public API
 * ----------
 *  GF.*                         — field arithmetic helpers
 *  generatorPoly(nsym)          — build generator polynomial for nsym ECC symbols
 *  berlekampMassey(syndromes)   — error-locator polynomial
 *  chienSearch(locator, n)      — error positions
 *  forney(syndromes, locator, positions, n) — error magnitudes
 *  erasureLocator(positions, n) — locator for known-erasure positions
 *  rsEncode(msg, nsym)          — low-level encode → Uint8Array
 *  rsDecode(msg, nsym, erasures) — low-level decode → { data, errata, corrections }
 *  RSCodec                      — high-level class
 */

"use strict";

// ─────────────────────────────────────────────────────────────────────────────
// GF(256) field arithmetic
// ─────────────────────────────────────────────────────────────────────────────

const GF_EXP = new Uint8Array(512);   // anti-log (exp) table, doubled to avoid modulo
const GF_LOG = new Uint8Array(256);   // log table

(function buildTables() {
  let x = 1;
  for (let i = 0; i < 255; i++) {
    GF_EXP[i] = x;
    GF_LOG[x] = i;
    x <<= 1;
    if (x & 0x100) x ^= 0x11d;       // reduce by primitive polynomial
  }
  for (let i = 255; i < 512; i++) {
    GF_EXP[i] = GF_EXP[i - 255];
  }
})();

/** GF(256) field helpers */
export const GF = {
  /** Addition (= subtraction) in GF(256) */
  add(a, b) { return a ^ b; },

  /** Subtraction — identical to addition in characteristic-2 fields */
  sub(a, b) { return a ^ b; },

  /** Multiplication in GF(256) */
  mul(a, b) {
    if (a === 0 || b === 0) return 0;
    return GF_EXP[(GF_LOG[a] + GF_LOG[b]) % 255];
  },

  /** Division in GF(256).  Throws on division by zero. */
  div(a, b) {
    if (b === 0) throw new RangeError("GF division by zero");
    if (a === 0) return 0;
    return GF_EXP[(GF_LOG[a] - GF_LOG[b] + 255) % 255];
  },

  /** Multiplicative inverse.  Throws for 0. */
  inverse(a) {
    if (a === 0) throw new RangeError("GF inverse of zero");
    return GF_EXP[255 - GF_LOG[a]];
  },

  /** a raised to the power n */
  pow(a, n) {
    return GF_EXP[(GF_LOG[a] * n) % 255];
  },

  // ── Polynomial helpers (polynomials are plain arrays, index 0 = highest degree) ──

  /** Scale every coefficient of poly by scalar */
  polyScale(poly, x) {
    return poly.map(c => GF.mul(c, x));
  },

  /** Add (XOR) two polynomials */
  polyAdd(p, q) {
    const out = new Array(Math.max(p.length, q.length)).fill(0);
    for (let i = 0; i < p.length; i++) out[i + out.length - p.length] ^= p[i];
    for (let i = 0; i < q.length; i++) out[i + out.length - q.length] ^= q[i];
    return out;
  },

  /** Multiply two polynomials */
  polyMul(p, q) {
    const out = new Array(p.length + q.length - 1).fill(0);
    for (let j = 0; j < q.length; j++) {
      for (let i = 0; i < p.length; i++) {
        out[i + j] ^= GF.mul(p[i], q[j]);
      }
    }
    return out;
  },

  /** Divide poly p by poly q → { quotient, remainder } */
  polyDiv(p, q) {
    let rem = [...p];
    const lead = q[0];
    const qdeg = q.length - 1;
    for (let i = 0; i < p.length - qdeg; i++) {
      const coef = GF.div(rem[i], lead);
      if (coef !== 0) {
        for (let j = 1; j < q.length; j++) {
          rem[i + j] ^= GF.mul(coef, q[j]);
        }
      }
    }
    const sep = p.length - qdeg;
    return {
      quotient: rem.slice(0, sep),
      remainder: rem.slice(sep),
    };
  },

  /** Evaluate polynomial at x using Horner's method */
  polyEval(poly, x) {
    let y = poly[0];
    for (let i = 1; i < poly.length; i++) {
      y = GF.add(GF.mul(y, x), poly[i]);
    }
    return y;
  },
};

// ─────────────────────────────────────────────────────────────────────────────
// Generator polynomial
// ─────────────────────────────────────────────────────────────────────────────

/**
 * Build the generator polynomial for a RS code with `nsym` ECC symbols.
 *
 * g(x) = ∏_{i=0}^{nsym-1} (x − α^i)
 *
 * @param {number} nsym  Number of ECC (parity) symbols.
 * @returns {number[]}   Coefficients, highest-degree first.
 */
export function generatorPoly(nsym) {
  let g = [1];
  for (let i = 0; i < nsym; i++) {
    g = GF.polyMul(g, [1, GF_EXP[i]]);
  }
  return g;
}

// ─────────────────────────────────────────────────────────────────────────────
// Syndrome computation  (internal helper)
// ─────────────────────────────────────────────────────────────────────────────

function calcSyndromes(msg, nsym) {
  const syn = [];
  for (let i = 0; i < nsym; i++) {
    syn.push(GF.polyEval(msg, GF_EXP[i]));
  }
  return syn;
}

// ─────────────────────────────────────────────────────────────────────────────
// Berlekamp-Massey  — error-locator polynomial
// ─────────────────────────────────────────────────────────────────────────────

/**
 * Berlekamp-Massey algorithm.
 *
 * Given the syndrome vector, find the shortest LFSR (error-locator polynomial σ)
 * that generates those syndromes.
 *
 * @param {number[]} syndromes
 * @returns {number[]} Error-locator polynomial (highest degree first).
 */
export function berlekampMassey(syndromes) {
  let C = [1];   // current connection polynomial
  let B = [1];   // previous connection polynomial
  let L = 0;
  let m = 1;
  let b = 1;

  for (let n = 0; n < syndromes.length; n++) {
    // Discrepancy
    let d = syndromes[n];
    for (let i = 1; i <= L; i++) {
      d ^= GF.mul(C[C.length - 1 - i] ?? 0, syndromes[n - i]);
    }

    if (d === 0) {
      m++;
    } else if (2 * L <= n) {
      const T = [...C];
      const scale = GF.div(d, b);
      // C = C − (d/b) · x^m · B
      const xmB = new Array(m).fill(0).concat(B);
      const scaled = GF.polyScale(xmB, scale);
      C = GF.polyAdd(C, scaled);
      L = n + 1 - L;
      B = T;
      b = d;
      m = 1;
    } else {
      const scale = GF.div(d, b);
      const xmB = new Array(m).fill(0).concat(B);
      C = GF.polyAdd(C, GF.polyScale(xmB, scale));
      m++;
    }
  }
  return C;
}

// ─────────────────────────────────────────────────────────────────────────────
// Chien search  — error positions
// ─────────────────────────────────────────────────────────────────────────────

/**
 * Chien search: find all roots of the error-locator polynomial in GF(256).
 *
 * A root α^{-i} means there is an error at position (n − 1 − i) counting from
 * the start of the codeword of length `n`.
 *
 * @param {number[]} locator  Error-locator polynomial (highest-degree first).
 * @param {number}   n        Codeword length.
 * @returns {number[]}        Error positions (byte indices within the codeword).
 */
export function chienSearch(locator, n) {
  const positions = [];
  for (let i = 0; i < n; i++) {
    if (GF.polyEval(locator, GF_EXP[255 - i]) === 0) {
      positions.push(n - 1 - i);
    }
  }
  return positions;
}

// ─────────────────────────────────────────────────────────────────────────────
// Forney algorithm  — error magnitudes
// ─────────────────────────────────────────────────────────────────────────────

/**
 * Forney algorithm: compute error magnitudes at known positions.
 *
 * @param {number[]} syndromes   Raw syndromes array (length nsym).
 * @param {number[]} locator     Error-locator polynomial σ(x).
 * @param {number[]} positions   Error positions (byte indices).
 * @param {number}   n           Codeword length.
 * @returns {number[]}           Error magnitudes, parallel to `positions`.
 */
export function forney(syndromes, locator, positions, n) {
  // Error evaluator polynomial Ω = (S · σ) mod x^nsym
  const nsym = syndromes.length;
  const synPoly = [1, ...syndromes];          // prepend leading 1 (degree-0 term)
  const raw = GF.polyMul(synPoly, locator);
  const omega = raw.slice(raw.length - nsym); // keep only last nsym coefficients

  // Formal derivative of σ (in char-2: drop even-index terms)
  const sigmaPrime = locator.slice().reverse().reduce((acc, c, i) => {
    if (i % 2 === 1) acc.push(c);
    return acc;
  }, []).reverse();

  const magnitudes = positions.map(pos => {
    const xi = GF_EXP[pos];                          // α^pos
    const xiInv = GF.inverse(xi);
    const num = GF.mul(xi, GF.polyEval(omega, xiInv));
    const den = GF.polyEval(sigmaPrime.length ? sigmaPrime : [1], xiInv);
    return GF.div(num, den);
  });

  return magnitudes;
}

// ─────────────────────────────────────────────────────────────────────────────
// Erasure locator  — for known-position damage (e.g. strand breaks)
// ─────────────────────────────────────────────────────────────────────────────

/**
 * Build the erasure-locator polynomial Γ(x) = ∏ (1 − α^{pos_i} · x).
 *
 * When erasure positions are known in advance (e.g. a broken DNA strand
 * segment whose index is known), this polynomial seeds the BM step and
 * lets you correct up to nsym erasures instead of only nsym/2 errors.
 *
 * @param {number[]} positions  Erasure positions (byte indices, 0-based).
 * @param {number}   n          Codeword length (unused here; kept for API symmetry).
 * @returns {number[]}          Erasure-locator polynomial (highest-degree first).
 */
export function erasureLocator(positions, n) {
  let elp = [1];
  for (const pos of positions) {
    elp = GF.polyMul(elp, GF.polyAdd([1], [GF_EXP[pos], 0]));
  }
  return elp;
}

// ─────────────────────────────────────────────────────────────────────────────
// rsEncode  — low-level encode
// ─────────────────────────────────────────────────────────────────────────────

/**
 * Encode a message with Reed-Solomon ECC.
 *
 * The ECC symbols are *appended* to the message (systematic form):
 *   codeword = [ msg... , ecc... ]
 *
 * @param {Uint8Array|number[]} msg   Raw data bytes.
 * @param {number}              nsym  Number of ECC symbols to append.
 * @returns {Uint8Array}              Full codeword (msg.length + nsym bytes).
 */
export function rsEncode(msg, nsym) {
  const gen = generatorPoly(nsym);
  const msgOut = new Uint8Array(msg.length + nsym);
  msgOut.set(msg);

  for (let i = 0; i < msg.length; i++) {
    const coef = msgOut[i];
    if (coef !== 0) {
      for (let j = 1; j < gen.length; j++) {
        msgOut[i + j] ^= GF.mul(gen[j], coef);
      }
    }
  }

  // Copy original message back (polynomial division corrupts leading bytes)
  const out = new Uint8Array(msg.length + nsym);
  out.set(msg);
  out.set(msgOut.slice(msg.length), msg.length);
  return out;
}

// ─────────────────────────────────────────────────────────────────────────────
// rsDecode  — low-level decode
// ─────────────────────────────────────────────────────────────────────────────

/**
 * Decode (and correct) a possibly-corrupted codeword.
 *
 * @param {Uint8Array|number[]} msg       Received codeword (data + ECC).
 * @param {number}              nsym      Number of ECC symbols.
 * @param {number[]}            [erasures=[]]  Known erasure positions (byte indices).
 * @returns {{ data: Uint8Array, errata: number[], corrections: number }}
 *   - data        Corrected data bytes (without ECC).
 *   - errata      Positions that were corrected.
 *   - corrections Number of symbol corrections applied.
 * @throws {Error} If the message is uncorrectable.
 */
export function rsDecode(msg, nsym, erasures = []) {
  if (erasures.length > nsym) {
    throw new Error("Too many erasures: " + erasures.length + " > nsym " + nsym);
  }

  const msgArr = Array.from(msg);
  const n = msgArr.length;

  // ── 1. Compute syndromes ──
  const syndromes = calcSyndromes(msgArr, nsym);
  const allZero = syndromes.every(s => s === 0);
  if (allZero) {
    return {
      data: new Uint8Array(msgArr.slice(0, n - nsym)),
      errata: [],
      corrections: 0,
    };
  }

  // ── 2. Erasure locator (if any erasures provided) ──
  const elp = erasures.length ? erasureLocator(erasures, n) : [1];

  // ── 3. Modified syndromes for combined erasure+error correction ──
  // S'(x) = S(x) · Γ(x)  mod x^nsym
  let synShifted = [...syndromes];
  for (const epos of erasures) {
    const xi = GF_EXP[epos];
    for (let j = 0; j < nsym - 1; j++) {
      synShifted[j] = GF.add(GF.mul(synShifted[j], xi), synShifted[j + 1]);
    }
    synShifted[nsym - 1] = GF.mul(synShifted[nsym - 1], xi);
  }

  // ── 4. BM on (possibly modified) syndromes for additional errors ──
  const errLocator = berlekampMassey(synShifted.slice(0, nsym - erasures.length));
  const numErrors = errLocator.length - 1;

  if (2 * numErrors + erasures.length > nsym) {
    throw new Error(
      `Message is uncorrectable: ${numErrors} errors + ${erasures.length} erasures > nsym ${nsym}`
    );
  }

  // ── 5. Full errata locator = erasure locator × error locator ──
  const errataLocator = GF.polyMul(elp, errLocator);

  // ── 6. Find all errata positions via Chien search ──
  const errataPositions = chienSearch(errataLocator, n);

  if (errataPositions.length !== errataLocator.length - 1) {
    throw new Error("Chien search found incorrect number of roots; message may be uncorrectable.");
  }

  // ── 7. Forney — compute magnitudes ──
  const magnitudes = forney(syndromes, errataLocator, errataPositions, n);

  // ── 8. Apply corrections ──
  const corrected = [...msgArr];
  for (let i = 0; i < errataPositions.length; i++) {
    corrected[errataPositions[i]] ^= magnitudes[i];
  }

  // ── 9. Verify ──
  const checkSyn = calcSyndromes(corrected, nsym);
  if (!checkSyn.every(s => s === 0)) {
    throw new Error("Decoding failed: residual syndromes are non-zero after correction.");
  }

  return {
    data: new Uint8Array(corrected.slice(0, n - nsym)),
    errata: errataPositions,
    corrections: errataPositions.length,
  };
}

// ─────────────────────────────────────────────────────────────────────────────
// RSCodec  — high-level class
// ─────────────────────────────────────────────────────────────────────────────

/**
 * High-level Reed-Solomon codec.
 *
 * @example
 * const rs = new RSCodec(10);              // 10 ECC symbols
 * const codeword = rs.encode("hello");     // Uint8Array
 * const { data } = rs.decode(codeword);   // original bytes
 *
 * // DNA strand workflow
 * const strand = rs.encodeStrand(data);    // base-4 string "ACGT…"
 * const decoded = rs.decodeStrand(strand); // { data, errata, corrections }
 */
export class RSCodec {
  /**
   * @param {number} nsym   ECC symbol count (default 10).
   *                        Corrects up to ⌊nsym/2⌋ errors, or up to nsym erasures.
   * @param {number} [blockSize=255]  Max codeword length (≤ 255 for GF(256)).
   */
  constructor(nsym = 10, blockSize = 255) {
    if (nsym < 1 || nsym > 254) throw new RangeError("nsym must be 1–254");
    if (blockSize < nsym + 1 || blockSize > 255) {
      throw new RangeError("blockSize must be nsym+1 to 255");
    }
    this.nsym = nsym;
    this.blockSize = blockSize;
    this._dataSize = blockSize - nsym;
    this._gen = generatorPoly(nsym);
    Object.freeze(this);
  }

  // ── Capacity ───────────────────────────────────────────────────────────────

  /**
   * Maximum correctable errors and erasures for this codec configuration.
   *
   * @returns {{ maxErrors: number, maxErasures: number, dataPerBlock: number, eccPerBlock: number }}
   */
  capacity() {
    return {
      maxErrors:    Math.floor(this.nsym / 2),
      maxErasures:  this.nsym,
      dataPerBlock: this._dataSize,
      eccPerBlock:  this.nsym,
      blockSize:    this.blockSize,
    };
  }

  // ── Core encode / decode ───────────────────────────────────────────────────

  /**
   * Encode data, splitting into blocks if necessary.
   *
   * @param {string|Uint8Array|number[]} data
   * @returns {Uint8Array}  Concatenated encoded blocks.
   */
  encode(data) {
    const bytes = this._toBytes(data);
    const chunks = this._chunkify(bytes, this._dataSize);
    const parts = chunks.map(chunk => rsEncode(chunk, this.nsym));
    return this._concat(parts);
  }

  /**
   * Decode and correct a received codeword (or multi-block sequence).
   *
   * @param {string|Uint8Array|number[]} data        Received (possibly corrupted) codeword.
   * @param {number[][]}                [erasures=[]] Per-block erasure positions.
   *   Pass a flat array for a single block, or a nested array for multi-block.
   * @returns {{ data: Uint8Array, errata: number[][], corrections: number }}
   */
  decode(data, erasures = []) {
    const bytes = this._toBytes(data);
    const chunks = this._chunkify(bytes, this.blockSize);

    const normalizeErasures = (e, blockIdx) => {
      if (!e.length) return [];
      if (Array.isArray(e[0])) return e[blockIdx] ?? [];
      return blockIdx === 0 ? e : [];
    };

    const decoded = chunks.map((chunk, i) =>
      rsDecode(chunk, this.nsym, normalizeErasures(erasures, i))
    );

    return {
      data:        this._concat(decoded.map(d => d.data)),
      errata:      decoded.map(d => d.errata),
      corrections: decoded.reduce((s, d) => s + d.corrections, 0),
    };
  }

  // ── Strand (DNA / base-4) codec ────────────────────────────────────────────

  /**
   * Encode data into a DNA-style base-4 string ("A", "C", "G", "T").
   *
   * Each byte is encoded as two nucleotides (high nibble, low nibble).
   *
   * @param {string|Uint8Array|number[]} data
   * @returns {string}  Nucleotide string of length (encoded bytes × 2).
   */
  encodeStrand(data) {
    const codeword = this.encode(data);
    return this._bytesToStrand(codeword);
  }

  /**
   * Decode a nucleotide strand back to the original data.
   *
   * @param {string}    strand              DNA nucleotide string.
   * @param {number[]}  [erasedNucleotides=[]] Positions of damaged/unknown nucleotides.
   *   Nucleotide positions are converted to byte positions automatically.
   * @returns {{ data: Uint8Array, errata: number[][], corrections: number }}
   */
  decodeStrand(strand, erasedNucleotides = []) {
    const bytes = this._strandToBytes(strand);
    const erasedBytes = this._nucleotidePosToBytePos(erasedNucleotides);
    return this.decode(bytes, erasedBytes);
  }

  // ── Internal helpers ───────────────────────────────────────────────────────

  /** Accept string, Uint8Array, or plain Array; always return Uint8Array. */
  _toBytes(data) {
    if (typeof data === "string") {
      return new TextEncoder().encode(data);
    }
    return data instanceof Uint8Array ? data : new Uint8Array(data);
  }

  /** Split a Uint8Array into chunks of `size`. */
  _chunkify(data, size) {
    const chunks = [];
    for (let i = 0; i < data.length; i += size) {
      chunks.push(data.slice(i, i + size));
    }
    return chunks.length ? chunks : [new Uint8Array(0)];
  }

  /** Concatenate an array of Uint8Arrays. */
  _concat(arrays) {
    const total = arrays.reduce((s, a) => s + a.length, 0);
    const out = new Uint8Array(total);
    let offset = 0;
    for (const arr of arrays) {
      out.set(arr, offset);
      offset += arr.length;
    }
    return out;
  }

  /** Convert bytes to a base-4 nucleotide string (A/C/G/T). */
  _bytesToStrand(bytes) {
    const BASES = "ACGT";
    let s = "";
    for (const byte of bytes) {
      s += BASES[(byte >> 6) & 0x3];
      s += BASES[(byte >> 4) & 0x3];
      s += BASES[(byte >> 2) & 0x3];
      s += BASES[byte & 0x3];
    }
    return s;
  }

  /** Convert a nucleotide string back to bytes; unknown chars become 0. */
  _strandToBytes(strand) {
    const MAP = { A: 0, C: 1, G: 2, T: 3 };
    const len = Math.ceil(strand.length / 4);
    const out = new Uint8Array(len);
    for (let i = 0; i < len; i++) {
      const b0 = MAP[strand[i * 4]]     ?? 0;
      const b1 = MAP[strand[i * 4 + 1]] ?? 0;
      const b2 = MAP[strand[i * 4 + 2]] ?? 0;
      const b3 = MAP[strand[i * 4 + 3]] ?? 0;
      out[i] = (b0 << 6) | (b1 << 4) | (b2 << 2) | b3;
    }
    return out;
  }

  /**
   * Map nucleotide (base-4 char) positions to byte positions.
   * Every 4 nucleotides correspond to 1 byte.
   * Multiple erased nucleotides in the same byte mark that byte erased once.
   */
  _nucleotidePosToBytePos(nucleotidePositions) {
    const byteSet = new Set(nucleotidePositions.map(p => Math.floor(p / 4)));
    return [...byteSet].sort((a, b) => a - b);
  }
}

// ─────────────────────────────────────────────────────────────────────────────
// Default export
// ─────────────────────────────────────────────────────────────────────────────

export default RSCodec;