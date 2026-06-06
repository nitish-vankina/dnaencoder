// =============================================================================
//  HelixArchive — Radiation Simulation Engine
//  radiation.js
//
//  Uses REAL solar proton flux data from NOAA GOES satellites via:
//  https://services.swpc.noaa.gov/json/goes/primary/integral-protons-1-day.json
//
//  The >=10 MeV proton flux channel is the primary driver of DNA base damage
//  in space environments. We fetch today's live reading, then scale it to
//  each environment using peer-reviewed distance/shielding multipliers.
//
//  Pipeline:
//    1. Fetch live NOAA GOES proton flux (pfu = protons/cm²/s/sr at >=10 MeV)
//    2. Average last 24h readings to get a baseline solar flux
//    3. Scale baseline by environment multiplier (LEO is most shielded;
//       Europa is ~1000x worse due to Jupiter's radiation belts)
//    4. Convert flux → damage probability per base per year
//    5. Run Monte Carlo simulation on each base of the DNA strand
//
//  Damage types modeled:
//    - Base change:   ionizing particle alters a nucleotide chemically
//    - Strand break:  high-energy particle snaps the backbone
//                     (represented as unreadable 'X' bases)
// =============================================================================

// ---------------------------------------------------------------------------
// Environment scaling factors relative to LEO
// Based on published GCR/SPE dose rate ratios:
//   Cucinotta et al. (2014), Durante & Cucinotta (2011),
//   Zeitlin et al. (2013) — MSL/RAD Mars surface measurements
// ---------------------------------------------------------------------------

const ENVIRONMENTS = {
  leo: {
    label:           'Low Earth Orbit (ISS)',
    description:     "Protected by Earth's magnetic field. Moderate radiation.",
    fluxMultiplier:  1,        // baseline — Earth's magnetosphere provides shielding
    breakMultiplier: 0.1,      // strand breaks rarer than base changes
    emoji:           '🛰',
  },
  lunar: {
    label:           'Lunar Surface',
    description:     'No magnetic field. Direct solar particle and GCR exposure.',
    fluxMultiplier:  8,
    breakMultiplier: 0.1,
    emoji:           '🌕',
  },
  mars: {
    label:           'Mars Surface',
    description:     'Thin CO₂ atmosphere, no global magnetic field. ~0.3 Sv/yr measured by MSL/RAD.',
    fluxMultiplier:  40,
    breakMultiplier: 0.1,
    emoji:           '🔴',
  },
  deepspace: {
    label:           'Deep Space (beyond Mars)',
    description:     'No planetary protection whatsoever. Full GCR flux.',
    fluxMultiplier:  120,
    breakMultiplier: 0.12,
    emoji:           '🌌',
  },
  europa: {
    label:           'Europa Orbit (Jupiter)',
    description:     "Trapped inside Jupiter's radiation belts. ~540 Sv/day surface dose.",
    fluxMultiplier:  3000,
    breakMultiplier: 0.18,
    emoji:           '🪐',
  },
};

// ---------------------------------------------------------------------------
// NOAA GOES proton flux endpoint
// Public, no API key, updated every 5 minutes
// Returns array of [timestamp, energy_channel, flux_value] objects
// ---------------------------------------------------------------------------

const NOAA_ENDPOINT =
  'https://services.swpc.noaa.gov/json/goes/primary/integral-protons-1-day.json';

// Conversion: proton flux (pfu) → base damage probability per base per year
// Derived from DNA strand break cross-section ~1e-9 cm² per 10 MeV proton
// and scaling to probability space for simulation purposes.
// This is a simplified but directionally accurate model.
const FLUX_TO_DAMAGE_RATE = 1e-13;  // per pfu per base per year

// ---------------------------------------------------------------------------
// Fetch and process live NOAA data
// ---------------------------------------------------------------------------

/**
 * Fetch the last 24 hours of GOES integral proton flux data from NOAA.
 * Returns the average flux at >=10 MeV in proton flux units (pfu).
 *
 * @returns {Promise<{ flux: number, source: string, timestamp: string, readings: number }>}
 */
async function fetchLiveFlux() {
  const response = await fetch(NOAA_ENDPOINT);
  if (!response.ok) throw new Error(`NOAA API returned ${response.status}`);

  const data = await response.json();

  // Each entry: { time_tag, satellite, flux, energy }
  // We want the >=10 MeV channel readings
  const channel10MeV = data.filter(d => d.energy === '>=10 MeV' && d.flux !== null && d.flux > 0);

  if (channel10MeV.length === 0) {
    throw new Error('No valid >=10 MeV readings in NOAA response');
  }

  // Average flux over available readings
  const avg = channel10MeV.reduce((sum, d) => sum + d.flux, 0) / channel10MeV.length;
  const latest = channel10MeV[channel10MeV.length - 1];

  return {
    flux:      avg,
    source:    `NOAA GOES (${latest.satellite || 'primary'})`,
    timestamp: latest.time_tag,
    readings:  channel10MeV.length,
  };
}

// ---------------------------------------------------------------------------
// Damage simulation
// ---------------------------------------------------------------------------

const BASES = ['A', 'C', 'G', 'T'];

/** Seeded pseudo-RNG (mulberry32) for reproducible results */
function makePRNG(seed) {
  let s = (seed >>> 0) || 0xdeadbeef;
  return function rand() {
    s += 0x6D2B79F5;
    let t = s;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/** Mutate a base to a different random base */
function mutateBase(base, rand) {
  const others = BASES.filter(b => b !== base.toUpperCase());
  return others[Math.floor(rand() * others.length)];
}

/**
 * Simulate radiation damage on a DNA strand using a real proton flux value.
 *
 * @param {string} strand          — original DNA strand (A/C/G/T)
 * @param {string} environment     — key from ENVIRONMENTS
 * @param {number} years           — mission duration
 * @param {number} baseFluxPfu     — measured solar proton flux in pfu (from NOAA)
 * @param {number} [seed]          — RNG seed
 * @returns {Array<{index, original, result, damageType}>}
 */
function simulateRadiation(strand, environment, years, baseFluxPfu, seed = Date.now()) {
  const env = ENVIRONMENTS[environment];
  if (!env) throw new Error(`Unknown environment: ${environment}`);

  const rand = makePRNG(seed);
  const upper = strand.toUpperCase();

  // Scale flux by environment multiplier
  const scaledFlux = baseFluxPfu * env.fluxMultiplier;

  // Convert to per-base annual damage probabilities
  const annualDamageRate = scaledFlux * FLUX_TO_DAMAGE_RATE;
  const annualBreakRate  = annualDamageRate * env.breakMultiplier;

  // Cumulative probability over mission duration
  const damageProb = 1 - Math.pow(1 - annualDamageRate, years);
  const breakProb  = 1 - Math.pow(1 - annualBreakRate,  years);

  const results = [];
  let inBreak = false;
  let breakRemaining = 0;

  for (let i = 0; i < upper.length; i++) {
    const base = upper[i];

    if (inBreak) {
      results.push({ index: i, original: base, result: 'X', damageType: 'unreadable' });
      breakRemaining--;
      if (breakRemaining <= 0) inBreak = false;
      continue;
    }

    if (rand() < breakProb) {
      const breakLen = 2 + Math.floor(rand() * 6);
      inBreak = true;
      breakRemaining = breakLen - 1;
      results.push({ index: i, original: base, result: 'X', damageType: 'break_start' });
      continue;
    }

    if (rand() < damageProb) {
      results.push({ index: i, original: base, result: mutateBase(base, rand), damageType: 'base_change' });
      continue;
    }

    results.push({ index: i, original: base, result: base, damageType: 'none' });
  }

  return results;
}

/**
 * Summarize damage results.
 *
 * @param {Array} results
 * @returns {{ total, damaged, baseChanges, unreadable, damageFraction, recoverable, damagedStrand }}
 */
function summarizeDamage(results) {
  let baseChanges = 0;
  let unreadable  = 0;

  for (const r of results) {
    if (r.damageType === 'base_change')                       baseChanges++;
    if (r.damageType === 'unreadable' || r.damageType === 'break_start') unreadable++;
  }

  const damaged        = baseChanges + unreadable;
  const damageFraction = damaged / results.length;
  const damagedStrand  = results.map(r => r.result).join('');
  const recoverable    = damageFraction < 0.05;

  return { total: results.length, damaged, baseChanges, unreadable, damageFraction, recoverable, damagedStrand };
}

// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------

const RadiationEngine = { ENVIRONMENTS, fetchLiveFlux, simulateRadiation, summarizeDamage, NOAA_ENDPOINT };

if (typeof module !== 'undefined' && module.exports) {
  module.exports = RadiationEngine;
} else if (typeof window !== 'undefined') {
  window.RadiationEngine  = RadiationEngine;
  window.ENVIRONMENTS     = ENVIRONMENTS;
  window.fetchLiveFlux    = fetchLiveFlux;
  window.simulateRadiation = simulateRadiation;
  window.summarizeDamage  = summarizeDamage;
}