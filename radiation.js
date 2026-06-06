// =============================================================================
//  HelixArchive — Radiation Simulation Engine
//  radiation.js
//
//  Models the effect of ionizing space radiation on a DNA storage strand.
//
//  Real radiation effects on DNA:
//    - Base damage:   a base is chemically altered → treated as random substitution
//    - Strand breaks: backbone snaps → treated as a run of unreadable bases ('X')
//    - Cross-linking: bases bond incorrectly → treated as local substitution cluster
//
//  Radiation intensity is modeled on real cosmic ray / solar particle flux
//  at different orbital distances, expressed as events per base per year.
// =============================================================================


// ---------------------------------------------------------------------------
// Radiation environment presets
// (approximate ionizing event rates, normalized for modeling purposes)
// ---------------------------------------------------------------------------

const ENVIRONMENTS = {
  leo: {
    label: 'Low Earth Orbit (ISS)',
    description: 'Protected by Earth\'s magnetic field. Moderate radiation.',
    baseDamageRate:  0.00001,   // events per base per year
    strandBreakRate: 0.000001,
  },
  lunar: {
    label: 'Lunar Surface',
    description: 'No magnetic field. Higher solar particle exposure.',
    baseDamageRate:  0.00008,
    strandBreakRate: 0.000008,
  },
  mars: {
    label: 'Mars Surface',
    description: 'Thin atmosphere, no global magnetic field. Significant GCR exposure.',
    baseDamageRate:  0.0002,
    strandBreakRate: 0.00002,
  },
  deepspace: {
    label: 'Deep Space (beyond Mars)',
    description: 'No planetary protection. Full galactic cosmic ray exposure.',
    baseDamageRate:  0.0006,
    strandBreakRate: 0.00006,
  },
  europa: {
    label: 'Europa Orbit (Jupiter)',
    description: 'Trapped in Jupiter\'s massive radiation belts. Extreme environment.',
    baseDamageRate:  0.004,
    strandBreakRate: 0.0004,
  },
};

// ---------------------------------------------------------------------------
// Damage types
// ---------------------------------------------------------------------------

const DAMAGE = {
  NONE:        'none',
  BASE_CHANGE: 'base_change',   // base mutated to a different base
  UNREADABLE:  'unreadable',    // base destroyed — cannot be decoded
  BREAK_START: 'break_start',   // start of a strand break region
  BREAK_END:   'break_end',
};

const BASES = ['A', 'C', 'G', 'T'];

/**
 * Mutate a base to a different random base.
 * @param {string} base
 * @returns {string}
 */
function mutatBase(base) {
  const others = BASES.filter(b => b !== base.toUpperCase());
  return others[Math.floor(Math.random() * others.length)];
}


// ---------------------------------------------------------------------------
// Core simulation
// ---------------------------------------------------------------------------

/**
 * Simulate radiation damage on a DNA strand.
 *
 * Returns an array of per-base result objects, one per base in the strand.
 * Each object describes what happened to that base.
 *
 * @param {string} strand         — original DNA strand (A/C/G/T)
 * @param {string} environment    — key from ENVIRONMENTS
 * @param {number} years          — mission duration in years
 * @param {number} seed           — optional RNG seed for reproducibility
 * @returns {Array<{
 *   index:      number,
 *   original:   string,
 *   result:     string,
 *   damageType: string,
 * }>}
 */
function simulateRadiation(strand, environment, years, seed = Date.now()) {
  const env = ENVIRONMENTS[environment];
  if (!env) throw new Error(`Unknown environment: ${environment}`);

  // Simple seeded RNG (mulberry32)
  let s = seed >>> 0;
  function rand() {
    s += 0x6D2B79F5;
    let t = s;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  }

  const upper = strand.toUpperCase();
  const results = [];

  let inBreak = false;
  let breakRemaining = 0;

  for (let i = 0; i < upper.length; i++) {
    const base = upper[i];

    // Active strand break region
    if (inBreak) {
      results.push({ index: i, original: base, result: 'X', damageType: DAMAGE.UNREADABLE });
      breakRemaining--;
      if (breakRemaining <= 0) inBreak = false;
      continue;
    }

    // Check for strand break event
    const breakProb = 1 - Math.pow(1 - env.strandBreakRate, years);
    if (rand() < breakProb) {
      const breakLen = 2 + Math.floor(rand() * 6); // 2–7 bases destroyed
      inBreak = true;
      breakRemaining = breakLen - 1;
      results.push({ index: i, original: base, result: 'X', damageType: DAMAGE.BREAK_START });
      continue;
    }

    // Check for base damage event
    const damageProb = 1 - Math.pow(1 - env.baseDamageRate, years);
    if (rand() < damageProb) {
      const mutated = mutatBase(base);
      results.push({ index: i, original: base, result: mutated, damageType: DAMAGE.BASE_CHANGE });
      continue;
    }

    // No damage
    results.push({ index: i, original: base, result: base, damageType: DAMAGE.NONE });
  }

  return results;
}


// ---------------------------------------------------------------------------
// Summary statistics
// ---------------------------------------------------------------------------

/**
 * Summarize the damage results from simulateRadiation.
 *
 * @param {Array} results — output of simulateRadiation
 * @returns {{
 *   total:          number,
 *   damaged:        number,
 *   baseChanges:    number,
 *   unreadable:     number,
 *   damageFraction: number,
 *   recoverable:    boolean,
 *   damagedStrand:  string,
 * }}
 */
function summarizeDamage(results) {
  let baseChanges = 0;
  let unreadable  = 0;

  for (const r of results) {
    if (r.damageType === DAMAGE.BASE_CHANGE) baseChanges++;
    if (r.damageType === DAMAGE.UNREADABLE || r.damageType === DAMAGE.BREAK_START) unreadable++;
  }

  const damaged        = baseChanges + unreadable;
  const damageFraction = damaged / results.length;
  const damagedStrand  = results.map(r => r.result).join('');

  // Rough recoverability threshold — above ~5% damage, Reed-Solomon can't save it
  // without significant redundancy
  const recoverable = damageFraction < 0.05;

  return {
    total:          results.length,
    damaged,
    baseChanges,
    unreadable,
    damageFraction,
    recoverable,
    damagedStrand,
  };
}


// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------

const RadiationEngine = {
  ENVIRONMENTS,
  DAMAGE,
  simulateRadiation,
  summarizeDamage,
};

if (typeof module !== 'undefined' && module.exports) {
  module.exports = RadiationEngine;
} else if (typeof window !== 'undefined') {
  window.RadiationEngine = RadiationEngine;
  // Also expose globals used by inline scripts
  window.simulateRadiation = simulateRadiation;
  window.summarizeDamage   = summarizeDamage;
  window.ENVIRONMENTS      = ENVIRONMENTS;
  window.DAMAGE            = DAMAGE;
}
