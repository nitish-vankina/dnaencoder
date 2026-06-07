// =============================================================================
//  HelixArchive — GC Content Optimizer
//  gc-optimizer.js
//
//  Problem: Raw 2-bit encoding maps bytes directly to bases. Some text
//  produces strands with poor GC content — heavy in A/T or G/C — which
//  makes them harder to synthesize and store reliably.
//
//  Solution: XOR masking.
//
//  Before encoding, each byte is XOR'd with a mask byte. XOR is its own
//  inverse — applying the same mask again recovers the original. This means
//  the transformation is fully reversible without any extra metadata beyond
//  the mask itself.
//
//  How mask selection works:
//    1. Try all 256 possible mask values against the input bytes
//    2. For each mask, encode the XOR'd bytes and measure GC content
//    3. Pick the mask that brings GC closest to 50%
//    4. Store the mask as the first byte of the output strand (4 bases)
//    5. Decoder reads the mask from the first 4 bases, then reverses XOR
//
//  Result: the same data, encoded to a strand with near-ideal GC content,
//  with no loss of information and full reversibility.
//
//  Wire format:
//    [4 bases: mask byte] [4 bases per data byte ...]
//
//  Ideal GC target: 50% (± 10% acceptable range: 40–60%)
// =============================================================================

const OPT_GC_TARGET = 0.50;
const OPT_GC_MIN    = 0.40;
const OPT_GC_MAX    = 0.60;

const OPT_BASE_TABLE = ['A', 'C', 'G', 'T'];

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/**
 * Encode a single byte to 4 bases (same as engine.js encodeByte).
 * Duplicated here so gc-optimizer.js is self-contained.
 */
function _encodeByte(byte) {
  return (
    OPT_BASE_TABLE[(byte >> 6) & 0b11] +
    OPT_BASE_TABLE[(byte >> 4) & 0b11] +
    OPT_BASE_TABLE[(byte >> 2) & 0b11] +
    OPT_BASE_TABLE[ byte       & 0b11]
  );
}

/**
 * Compute GC fraction of a nucleotide string.
 */
function _gcFraction(strand) {
  if (!strand.length) return 0;
  let gc = 0;
  for (const ch of strand) {
    if (ch === 'G' || ch === 'C') gc++;
  }
  return gc / strand.length;
}

/**
 * Encode an array of bytes XOR'd with a given mask, return the strand.
 * Does not include the mask header — that's added by the caller.
 */
function _encodeWithMask(bytes, mask) {
  const out = [];
  for (let i = 0; i < bytes.length; i++) {
    out.push(_encodeByte(bytes[i] ^ mask));
  }
  return out.join('');
}


// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/**
 * Find the XOR mask (0–255) that brings GC content of the encoded strand
 * closest to the target (default 50%).
 *
 * Tries all 256 possible mask values and scores each by how close the
 * resulting GC is to OPT_GC_TARGET.
 *
 * @param   {Uint8Array} bytes
 * @returns {{ mask: number, gc: number, strand: string }}
 */
function findBestMask(bytes) {
  let bestMask  = 0;
  let bestGC    = _gcFraction(_encodeWithMask(bytes, 0));
  let bestDelta = Math.abs(bestGC - OPT_GC_TARGET);

  for (let mask = 1; mask <= 255; mask++) {
    const strand = _encodeWithMask(bytes, mask);
    const gc     = _gcFraction(strand);
    const delta  = Math.abs(gc - OPT_GC_TARGET);

    if (delta < bestDelta) {
      bestDelta = delta;
      bestMask  = mask;
      bestGC    = gc;
    }
  }

  return {
    mask:   bestMask,
    gc:     bestGC,
    strand: _encodeWithMask(bytes, bestMask),
  };
}

/**
 * Encode bytes to a GC-optimized DNA strand.
 *
 * Wire format: [4-base mask header] + [4 bases per data byte]
 * Total length: (bytes.length + 1) * 4
 *
 * @param   {Uint8Array} bytes   — raw UTF-8 bytes
 * @returns {{
 *   strand:      string,   — full optimized strand including mask header
 *   mask:        number,   — the XOR mask used (0–255)
 *   gcBefore:    number,   — GC fraction before optimization
 *   gcAfter:     number,   — GC fraction after optimization (data only)
 *   improvement: number,   — absolute GC improvement toward 50%
 *   inRange:     boolean,  — whether gcAfter is in [40%, 60%]
 * }}
 */
function optimizeGC(bytes) {
  // Measure unoptimized GC
  const rawStrand = _encodeWithMask(bytes, 0);
  const gcBefore  = _gcFraction(rawStrand);

  // Find the best mask
  const { mask, gc: gcAfter, strand: dataStrand } = findBestMask(bytes);

  // Prepend 4-base mask header
  const maskHeader = _encodeByte(mask);
  const fullStrand = maskHeader + dataStrand;

  const deltaBefore = Math.abs(gcBefore - OPT_GC_TARGET);
  const deltaAfter  = Math.abs(gcAfter  - OPT_GC_TARGET);

  return {
    strand:      fullStrand,
    mask,
    gcBefore,
    gcAfter,
    improvement: deltaBefore - deltaAfter,
    inRange:     gcAfter >= OPT_GC_MIN && gcAfter <= OPT_GC_MAX,
  };
}

/**
 * Decode a GC-optimized strand back to the original bytes.
 *
 * Reads the mask from the first 4 bases, then reverses the XOR on each
 * subsequent group of 4 bases.
 *
 * @param   {string} strand   — full optimized strand (including mask header)
 * @returns {Uint8Array}       original bytes
 * @throws  {Error}            if strand length is not a multiple of 4
 */
function decodeOptimized(strand) {
  if (strand.length % 4 !== 0) {
    throw new Error(
      `Strand length ${strand.length} is not a multiple of 4`
    );
  }

  const BASE_BITS = { A: 0b00, C: 0b01, G: 0b10, T: 0b11 };

  function decodeQuad(q) {
    let val = 0;
    for (let i = 0; i < 4; i++) {
      const bits = BASE_BITS[q[i].toUpperCase()];
      if (bits === undefined) throw new Error(`Invalid base '${q[i]}'`);
      val = (val << 2) | bits;
    }
    return val;
  }

  // First quad is the mask
  const mask = decodeQuad(strand.slice(0, 4));

  // Remaining quads are data bytes XOR'd with the mask
  const numBytes = (strand.length / 4) - 1;
  const out = new Uint8Array(numBytes);
  for (let i = 0; i < numBytes; i++) {
    out[i] = decodeQuad(strand.slice((i + 1) * 4, (i + 2) * 4)) ^ mask;
  }
  return out;
}

/**
 * Full pipeline: text → GC-optimized DNA strand.
 *
 * @param   {string} text
 * @returns {{ strand, mask, gcBefore, gcAfter, improvement, inRange }}
 */
function optimizeText(text) {
  const bytes = new TextEncoder().encode(text);
  return optimizeGC(bytes);
}

/**
 * Full pipeline: GC-optimized strand → original text.
 *
 * @param   {string} strand
 * @returns {string}
 */
function decodeOptimizedText(strand) {
  return new TextDecoder().decode(decodeOptimized(strand));
}


// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------

const GCOptimizer = {
  optimizeGC,
  optimizeText,
  decodeOptimized,
  decodeOptimizedText,
  findBestMask,
  OPT_GC_TARGET,
  OPT_GC_MIN,
  OPT_GC_MAX,
};

if (typeof module !== 'undefined' && module.exports) {
  module.exports = GCOptimizer;
} else if (typeof window !== 'undefined') {
  window.GCOptimizer          = GCOptimizer;
  window.optimizeGC           = optimizeGC;
  window.optimizeText         = optimizeText;
  window.decodeOptimized      = decodeOptimized;
  window.decodeOptimizedText  = decodeOptimizedText;
  window.findBestMask         = findBestMask;
}
