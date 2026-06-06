// =============================================================================
//  HelixArchive — Encoding Engine
//  engine.js
//
//  Pipeline:  Text  →  UTF-8 bytes  →  binary bits  →  DNA bases
//
//  Mapping (2 bits per base):
//    00 → A (Adenine)
//    01 → C (Cytosine)
//    10 → G (Guanine)
//    11 → T (Thymine)
//
//  Every byte produces 8 bits → 4 bases (4× expansion).
// =============================================================================


// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/** 2-bit → nucleotide lookup */
const BASE_TABLE = ['A', 'C', 'G', 'T'];

/** Nucleotide → 2-bit lookup */
const BASE_BITS = { A: 0b00, C: 0b01, G: 0b10, T: 0b11,
                    a: 0b00, c: 0b01, g: 0b10, t: 0b11 };

/** Max recommended consecutive identical bases for DNA synthesis */
const MAX_HOMOPOLYMER = 4;

/** Ideal GC content range for synthesis reliability */
const GC_MIN = 0.40;
const GC_MAX = 0.60;


// ---------------------------------------------------------------------------
// Step 1 — Text → UTF-8 bytes
// ---------------------------------------------------------------------------

/**
 * Encode a JavaScript string to a Uint8Array of UTF-8 bytes.
 *
 * @param  {string}     text
 * @returns {Uint8Array}
 *
 * @example
 * textToBytes("Hi")
 * // → Uint8Array [72, 105]
 */
function textToBytes(text) {
  return new TextEncoder().encode(text);
}

/**
 * Decode a UTF-8 Uint8Array back to a JavaScript string.
 *
 * @param  {Uint8Array} bytes
 * @returns {string}
 */
function bytesToText(bytes) {
  return new TextDecoder().decode(bytes);
}


// ---------------------------------------------------------------------------
// Step 2 — Bytes → binary bit string (diagnostic / display)
// ---------------------------------------------------------------------------

/**
 * Convert a Uint8Array to a binary string with each byte zero-padded to 8 bits.
 * Bytes are separated by a space for readability.
 *
 * This is purely for display — the encoding pipeline works directly on bytes,
 * not on a binary string, for performance.
 *
 * @param  {Uint8Array} bytes
 * @returns {string}
 *
 * @example
 * bytesToBinaryString(new Uint8Array([72]))
 * // → "01001000"
 */
function bytesToBinaryString(bytes) {
  return Array.from(bytes)
    .map(b => b.toString(2).padStart(8, '0'))
    .join(' ');
}

/**
 * Show how a single byte maps to its 4 diads and 4 bases.
 * Returns a structured object useful for step-by-step display.
 *
 * @param  {number} byte  — integer 0–255
 * @returns {{ byte, binary, diads, bases }}
 *
 * @example
 * explainByte(0x48)
 * // → { byte: 72, binary: "01001000",
 * //     diads: ["01","00","10","00"],
 * //     bases: ["C","A","G","A"] }
 */
function explainByte(byte) {
  const binary = byte.toString(2).padStart(8, '0');
  const diads  = [
    binary.slice(0, 2),
    binary.slice(2, 4),
    binary.slice(4, 6),
    binary.slice(6, 8),
  ];
  const bases = diads.map(d => BASE_TABLE[parseInt(d, 2)]);
  return { byte, binary, diads, bases };
}


// ---------------------------------------------------------------------------
// Step 3 — Bytes → DNA strand
// ---------------------------------------------------------------------------

/**
 * Encode a single byte (0–255) to exactly 4 nucleotide characters.
 *
 * The byte is split into four 2-bit diads, MSB first:
 *   byte = b7b6 | b5b4 | b3b2 | b1b0
 *
 * @param  {number} byte
 * @returns {string}  4-character nucleotide string
 *
 * @example
 * encodeByte(0x00) // → "AAAA"
 * encodeByte(0xFF) // → "TTTT"
 * encodeByte(0x48) // → "CAGA"  (H = 0x48 = 01 00 10 00)
 */
function encodeByte(byte) {
  return (
    BASE_TABLE[(byte >> 6) & 0b11] +
    BASE_TABLE[(byte >> 4) & 0b11] +
    BASE_TABLE[(byte >> 2) & 0b11] +
    BASE_TABLE[ byte       & 0b11]
  );
}

/**
 * Decode exactly 4 nucleotide characters back to a single byte.
 *
 * @param  {string} quad  — 4-character nucleotide string (case-insensitive)
 * @returns {number}       byte value 0–255
 * @throws  {Error}        if any character is not A/C/G/T
 *
 * @example
 * decodeQuad("CAGA") // → 72  (0x48 = 'H')
 */
function decodeQuad(quad) {
  let result = 0;
  for (let i = 0; i < 4; i++) {
    const bits = BASE_BITS[quad[i]];
    if (bits === undefined) {
      throw new Error(`Invalid nucleotide '${quad[i]}' at position ${i}`);
    }
    result = (result << 2) | bits;
  }
  return result;
}

/**
 * Encode a Uint8Array of bytes to a DNA strand string.
 * Each byte becomes 4 bases → output length = input.length × 4.
 *
 * @param  {Uint8Array} bytes
 * @returns {string}    nucleotide strand (uppercase A/C/G/T)
 *
 * @example
 * encodeBytes(new Uint8Array([72, 105]))
 * // → "CAGACGAC"
 */
function encodeBytes(bytes) {
  // Pre-allocate a char array for performance
  const out = new Array(bytes.length * 4);
  for (let i = 0; i < bytes.length; i++) {
    const b = bytes[i];
    out[i*4]   = BASE_TABLE[(b >> 6) & 0b11];
    out[i*4+1] = BASE_TABLE[(b >> 4) & 0b11];
    out[i*4+2] = BASE_TABLE[(b >> 2) & 0b11];
    out[i*4+3] = BASE_TABLE[ b       & 0b11];
  }
  return out.join('');
}

/**
 * Decode a DNA strand string back to a Uint8Array of bytes.
 * Strand length must be a multiple of 4.
 *
 * @param  {string}     strand  — nucleotide string (case-insensitive)
 * @returns {Uint8Array}
 * @throws  {Error}             on invalid bases or bad strand length
 *
 * @example
 * decodeStrand("CAGACGAC")
 * // → Uint8Array [72, 105]  →  "Hi"
 */
function decodeStrand(strand) {
  if (strand.length % 4 !== 0) {
    throw new Error(
      `Strand length ${strand.length} is not a multiple of 4 — cannot decode to whole bytes`
    );
  }
  const out = new Uint8Array(strand.length / 4);
  for (let i = 0; i < out.length; i++) {
    out[i] = decodeQuad(strand.slice(i * 4, i * 4 + 4));
  }
  return out;
}


// ---------------------------------------------------------------------------
// Full pipeline helpers
// ---------------------------------------------------------------------------

/**
 * Full encoding pipeline: text string → DNA strand string.
 *
 * Text → UTF-8 bytes → nucleotide bases
 *
 * @param  {string} text
 * @returns {string} DNA strand
 *
 * @example
 * encode("Hi") // → "CAGACGAC"
 */
function encode(text) {
  return encodeBytes(textToBytes(text));
}

/**
 * Full decoding pipeline: DNA strand string → original text string.
 *
 * Nucleotide bases → bytes → UTF-8 text
 *
 * @param  {string} strand
 * @returns {string} original text
 *
 * @example
 * decode("CAGACGAC") // → "Hi"
 */
function decode(strand) {
  return bytesToText(decodeStrand(strand));
}


// ---------------------------------------------------------------------------
// Biological constraint analysis
// ---------------------------------------------------------------------------

/**
 * Compute the GC content of a strand as a fraction [0, 1].
 *
 * @param  {string} strand
 * @returns {number}
 */
function gcContent(strand) {
  if (strand.length === 0) return 0;
  let gc = 0;
  for (const ch of strand) {
    if (ch === 'G' || ch === 'C' || ch === 'g' || ch === 'c') gc++;
  }
  return gc / strand.length;
}

/**
 * Find all homopolymer runs — consecutive identical bases — exceeding minLength.
 *
 * @param  {string} strand
 * @param  {number} minLength  — minimum run length to report (default: MAX_HOMOPOLYMER + 1 = 5)
 * @returns {Array<{ base, position, length }>}
 */
function findHomopolymerRuns(strand, minLength = MAX_HOMOPOLYMER + 1) {
  if (strand.length === 0) return [];
  const runs  = [];
  let start = 0;
  let len   = 1;
  const s   = strand.toUpperCase();

  for (let i = 1; i < s.length; i++) {
    if (s[i] === s[i - 1]) {
      len++;
    } else {
      if (len >= minLength) runs.push({ base: s[start], position: start, length: len });
      start = i;
      len   = 1;
    }
  }
  if (len >= minLength) runs.push({ base: s[start], position: start, length: len });
  return runs;
}

/**
 * Return the reverse complement of a DNA strand.
 *
 * A↔T, C↔G. Non-ACGT characters are preserved.
 *
 * @param  {string} strand
 * @returns {string}
 *
 * @example
 * reverseComplement("ACGT") // → "ACGT"
 * reverseComplement("AAAA") // → "TTTT"
 */
function reverseComplement(strand) {
  const comp = { A:'T', T:'A', C:'G', G:'C', a:'t', t:'a', c:'g', g:'c' };
  return strand.split('').reverse().map(ch => comp[ch] ?? ch).join('');
}

/**
 * Estimate melting temperature (Tm) in °C.
 *
 * Uses the SantaLucia 1998 nearest-neighbor thermodynamic model for strands
 * ≥ 14 bases (1 M NaCl, 250 nM strand concentration).
 * Falls back to the Wallace rule for shorter strands.
 *
 * @param  {string}      strand
 * @returns {number|null}        null if strand is empty or contains invalid bases
 */
function meltingTemp(strand) {
  if (strand.length === 0) return null;

  const upper = strand.toUpperCase();
  for (const ch of upper) {
    if (!'ACGT'.includes(ch)) return null;
  }

  // Wallace rule for short strands
  if (upper.length < 14) {
    const at = upper.split('').filter(b => b==='A'||b==='T').length;
    const gc = upper.split('').filter(b => b==='G'||b==='C').length;
    return 2 * at + 4 * gc;
  }

  // Nearest-neighbor parameters (SantaLucia 1998)
  // ΔH kcal/mol, ΔS cal/mol·K
  const NN = {
    AA:{h:-7.9,s:-22.2}, TT:{h:-7.9,s:-22.2},
    AT:{h:-7.2,s:-20.4}, TA:{h:-7.2,s:-21.3},
    CA:{h:-8.5,s:-22.7}, TG:{h:-8.5,s:-22.7},
    GT:{h:-8.4,s:-22.4}, AC:{h:-8.4,s:-22.4},
    CT:{h:-7.8,s:-21.0}, AG:{h:-7.8,s:-21.0},
    GA:{h:-8.2,s:-22.2}, TC:{h:-8.2,s:-22.2},
    CG:{h:-10.6,s:-27.2},GC:{h:-9.8,s:-24.4},
    GG:{h:-8.0,s:-19.9}, CC:{h:-8.0,s:-19.9},
  };

  let dh = 0, ds = 0;
  for (let i = 0; i < upper.length - 1; i++) {
    const pair = upper[i] + upper[i+1];
    if (NN[pair]) { dh += NN[pair].h; ds += NN[pair].s; }
  }

  // Initiation corrections
  const firstGC = upper[0] === 'G' || upper[0] === 'C';
  const lastGC  = upper[upper.length-1] === 'G' || upper[upper.length-1] === 'C';
  dh += firstGC ? 0.1 : 2.3;
  dh += lastGC  ? 0.1 : 2.3;
  ds += firstGC ? -2.8 : 4.1;
  ds += lastGC  ? -2.8 : 4.1;

  // Tm = ΔH / (ΔS + R·ln(CT/4)) − 273.15
  const R = 1.987, CT = 250e-9;
  return (dh * 1000) / (ds + R * Math.log(CT / 4)) - 273.15;
}

/**
 * Full strand analysis. Returns all biological metrics in one call.
 *
 * @param  {string} strand
 * @returns {{
 *   length:       number,
 *   gcFraction:   number,
 *   gcOk:         boolean,
 *   hpRuns:       Array,
 *   tm:           number|null,
 *   revComp:      string,
 *   baseCounts:   { A, C, G, T },
 * }}
 */
function analyzeStrand(strand) {
  const gc = gcContent(strand);
  const counts = { A:0, C:0, G:0, T:0 };
  for (const ch of strand.toUpperCase()) {
    if (counts[ch] !== undefined) counts[ch]++;
  }
  return {
    length:     strand.length,
    gcFraction: gc,
    gcOk:       gc >= GC_MIN && gc <= GC_MAX,
    hpRuns:     findHomopolymerRuns(strand),
    tm:         meltingTemp(strand.slice(0, 200)), // sample first 200 bases
    revComp:    reverseComplement(strand),
    baseCounts: counts,
  };
}


// ---------------------------------------------------------------------------
// Exports (works as ES module or plain <script> tag)
// ---------------------------------------------------------------------------

const HelixEngine = {
  // Pipeline
  encode,
  decode,
  // Step-by-step
  textToBytes,
  bytesToText,
  bytesToBinaryString,
  explainByte,
  encodeByte,
  decodeQuad,
  encodeBytes,
  decodeStrand,
  // Bio analysis
  gcContent,
  findHomopolymerRuns,
  reverseComplement,
  meltingTemp,
  analyzeStrand,
  // Constants
  BASE_TABLE,
  BASE_BITS,
  MAX_HOMOPOLYMER,
  GC_MIN,
  GC_MAX,
};

if (typeof module !== 'undefined' && module.exports) {
  module.exports = HelixEngine;                   // Node.js / CommonJS
} else if (typeof window !== 'undefined') {
  window.HelixEngine = HelixEngine;               // Browser <script> tag
}
// ---------------------------------------------------------------------------
// Global aliases (backwards compatibility)
// ---------------------------------------------------------------------------
// Some apps (like index.html) expect top-level functions. Provide safe
// aliases on `window` to avoid changing callers.
if (typeof window !== 'undefined') {
  window.textToBytes        = textToBytes;
  window.bytesToText        = bytesToText;
  window.encodeBytes        = encodeBytes;
  window.decodeStrand       = decodeStrand;
  window.gcContent          = gcContent;
  window.analyzeStrand      = analyzeStrand;
  window.encode             = encode;
}