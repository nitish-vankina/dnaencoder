const ENVIRONMENTS = {
  leo: {
    label: 'Low Earth Orbit (ISS)',
    description: 'Protected by Earth\'s magnetic field. Moderate radiation.',
    baseDamageRate:  0.00001,   // events per base per year
    strandBreakRate: 0.000001,
    fallbackFlux: 0.3,          // mSv/day — ISS average (~110 mSv/yr)
  },
  lunar: {
    label: 'Lunar Surface',
    description: 'No magnetic field. Higher solar particle exposure.',
    baseDamageRate:  0.00008,
    strandBreakRate: 0.000008,
    fallbackFlux: 1.37,         // mSv/day — LRO/CRaTER measurements (~0.5 Sv/yr surface)
  },
  mars: {
    label: 'Mars Surface',
    description: 'Thin atmosphere, no global magnetic field. Significant GCR exposure. ' +
                 'MSL/RAD instrument measured ~0.64 mSv/day on the surface (Hassler et al. 2014).',
    baseDamageRate:  0.0002,
    strandBreakRate: 0.00002,
    fallbackFlux: 0.64,         // mSv/day — MSL/RAD surface measurement (Hassler et al. 2014)
  },
  deepspace: {
    label: 'Deep Space (beyond Mars)',
    description: 'No planetary protection. Full galactic cosmic ray exposure.',
    baseDamageRate:  0.0006,
    strandBreakRate: 0.00006,
    fallbackFlux: 1.84,         // mSv/day — Curiosity cruise-phase RAD average (~0.67 Sv/yr)
  },
  europa: {
    label: 'Europa Orbit (Jupiter)',
    description: 'Trapped in Jupiter\'s massive radiation belts. Extreme environment.',
    baseDamageRate:  0.004,
    strandBreakRate: 0.0004,
    fallbackFlux: 5600,         // mSv/day — estimated Jovian belt flux at Europa's orbit
  },
};

// ---------------------------------------------------------------------------
// Live flux fetch
// ---------------------------------------------------------------------------

/**
 * Attempt to fetch a live radiation flux reading for the given environment.
 *
 * Currently only LEO has a public real-time data source (NASA's SpaceWeather
 * DONKI API).  All other environments fall back to the mission-dosimetry
 * constants stored in ENVIRONMENTS[env].fallbackFlux.
 *
 * @param {string} environment — key from ENVIRONMENTS
 * @returns {Promise<{ flux: number, source: 'live'|'fallback', units: string }>}
 */
async function fetchLiveFlux(environment) {
  const env = ENVIRONMENTS[environment];
  if (!env) throw new Error(`Unknown environment: ${environment}`);

  if (environment === 'leo') {
    try {
      // NASA DONKI solar energetic particle endpoint — no API key required for
      // basic queries.  We grab the most recent 7-day window and take the last
      // reported flux value.
      const today = new Date();
      const weekAgo = new Date(today.getTime() - 7 * 24 * 60 * 60 * 1000);
      const fmt = d => d.toISOString().slice(0, 10);
      const url =
        `https://kauai.ccmc.gsfc.nasa.gov/DONKI/WS/get/SEP` +
        `?startDate=${fmt(weekAgo)}&endDate=${fmt(today)}&speed=&halfAngle=&catalog=ALL&keyword=NONE`;

      const res = await fetch(url, { signal: AbortSignal.timeout(5000) });
      if (!res.ok) throw new Error(`DONKI HTTP ${res.status}`);

      const events = await res.json();
      if (Array.isArray(events) && events.length > 0) {
        // Pull the peak flux from the most recent event (units: pfu = particles/cm²/s/sr)
        const latest = events[events.length - 1];
        const flux = latest?.peakIntensity ?? null;
        if (flux !== null && isFinite(flux)) {
          return { flux, source: 'live', units: 'pfu' };
        }
      }
      // Fall through if no events in window
    } catch (_) {
      // Network error or timeout — fall through to fallback
    }
  }

  return { flux: env.fallbackFlux, source: 'fallback', units: 'mSv/day' };
}


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
function mutateBase(base) {
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
      const mutated = mutateBase(base);
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
  fetchLiveFlux,
  simulateRadiation,
  summarizeDamage,
};

if (typeof module !== 'undefined' && module.exports) {
  module.exports = RadiationEngine;
} else if (typeof window !== 'undefined') {
  window.RadiationEngine = RadiationEngine;
  // Also expose globals used by inline scripts
  window.fetchLiveFlux     = fetchLiveFlux;
  window.simulateRadiation = simulateRadiation;
  window.summarizeDamage   = summarizeDamage;
  window.ENVIRONMENTS      = ENVIRONMENTS;
  window.DAMAGE            = DAMAGE;
}
