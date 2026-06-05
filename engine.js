/*
 * ShockLogic engine — pure meta-analysis core for the IABP / cardiogenic-shock
 * mechanical-support dashboard. Hartung–Knapp–Sidik–Jonkman (HKSJ) pooling of
 * log odds ratios, plus a baseline-risk competing-risk projection.
 *
 * Extracted from the original inline web-worker so the statistical core is a
 * single source of truth, importable under Node for testing. Used by the page
 * synchronously (the worker was removed — the computation is trivially small).
 *
 * Two correctness fixes applied during the 2026-06 revival (see runHKSJ):
 *   1. HKSJ variance FLOOR — the original divided q/( (k-1)*wSum ) with no floor,
 *      so when Q < k-1 the HK interval narrowed *below* the random-effects
 *      interval (anti-conservative). Floored at the standard RE variance.
 *   2. Proper t critical value — the original hardcoded 2.78 (k<=5) / 2.26 else.
 *      Now uses qt(0.975, k-1): e.g. k=2 -> 12.71, k=3 -> 4.30, k=5 -> 2.78.
 */

// Inverse Student-t at the two-sided 0.975 level. Exact table for df 1..30
// (matches R qt() to 3 dp); Cornish–Fisher expansion for df > 30.
function tInv975(df) {
    const T = {
        1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365,
        8: 2.306, 9: 2.262, 10: 2.228, 11: 2.201, 12: 2.179, 13: 2.160, 14: 2.145,
        15: 2.131, 16: 2.120, 17: 2.110, 18: 2.101, 19: 2.093, 20: 2.086, 21: 2.080,
        22: 2.074, 23: 2.069, 24: 2.064, 25: 2.060, 26: 2.056, 27: 2.052, 28: 2.048,
        29: 2.045, 30: 2.042
    };
    if (df <= 0) return Infinity;
    const d = Math.round(df);
    if (d <= 30) return T[d];
    const z = 1.959963985;             // qnorm(0.975)
    const z2 = z * z, z3 = z2 * z, z5 = z3 * z2;
    return z
        + (z3 + z) / (4 * df)
        + (5 * z5 + 16 * z3 + 3 * z) / (96 * df * df);
}

// HKSJ pooling on the supplied effect/SE keys. studies: [{[kEff], [kSe], id}].
function runHKSJ(studies, kEff, kSe) {
    const k = studies.length;
    if (k === 1) {
        const s = studies[0];
        return {
            or: Math.exp(s[kEff]),
            lo: Math.exp(s[kEff] - 1.96 * s[kSe]),
            hi: Math.exp(s[kEff] + 1.96 * s[kSe]),
            I2: 0, tau2: 0, k: 1,
            studies: studies.map(st => oneStudy(st, kEff, kSe))
        };
    }

    const y = studies.map(s => s[kEff]);
    const v = studies.map(s => s[kSe] ** 2);

    // DerSimonian–Laird tau^2 (fixed-effect weights for Q and C)
    let num = 0, den = 0;
    v.forEach((val, i) => { num += (1 / val) * y[i]; den += (1 / val); });
    const muFE = num / den;

    let Q = 0;
    y.forEach((val, i) => Q += (1 / v[i]) * (val - muFE) ** 2);

    const C = den - (v.map(val => 1 / val ** 2).reduce((a, b) => a + b, 0) / den);
    const tau2 = Math.max(0, (Q - (k - 1)) / C);
    const I2 = Math.max(0, (Q - (k - 1)) / Q) * 100;

    // Random-effects weighted mean
    const wRE = v.map(val => 1 / (val + tau2));
    let wSum = 0, wMeanNum = 0;
    wRE.forEach((w, i) => { wSum += w; wMeanNum += w * y[i]; });
    const muRE = wMeanNum / wSum;

    // Hartung–Knapp variance with the standard floor at the RE variance.
    let qTerm = 0;
    wRE.forEach((w, i) => qTerm += w * (y[i] - muRE) ** 2);
    const qStat = qTerm / (k - 1);          // HK multiplier
    const varRE = 1 / wSum;                  // classical RE variance
    const varHKSJ = Math.max(1, qStat) * varRE;   // <-- floor (fix #1)
    const seHKSJ = Math.sqrt(varHKSJ);

    const tVal = tInv975(k - 1);             // <-- proper t_{k-1} (fix #2)

    return {
        or: Math.exp(muRE),
        lo: Math.exp(muRE - tVal * seHKSJ),
        hi: Math.exp(muRE + tVal * seHKSJ),
        I2, tau2, k,
        studies: studies.map(st => oneStudy(st, kEff, kSe))
    };
}

function oneStudy(st, kEff, kSe) {
    return {
        id: st.id,
        or: Math.exp(st[kEff]),
        lo: Math.exp(st[kEff] - 1.96 * st[kSe]),
        hi: Math.exp(st[kEff] + 1.96 * st[kSe])
    };
}

// Full simulation (ported verbatim from the old worker, minus message plumbing).
// Returns null when no device is active. devices: the DB rows with `active`.
function runShockSim({ devices, subgroup, timeDelay }) {
    const activeDevices = devices.filter(d => d.active);
    if (activeDevices.length === 0) return null;

    const adjustedData = activeDevices.map(d => {
        let logOR = Math.log(d.or);
        let logSafe = Math.log(d.or_safety);
        let se = d.se;

        if (subgroup === 'stemi') {
            if (d.id.includes('DanGer')) { logOR = Math.log(0.68); se = se * 0.9; }
            else if (d.id.includes('IABP')) { /* neutral */ }
            else { se = se * 1.2; }
        }

        if (logOR < 0) {
            const decay = timeDelay / 6;
            logOR = logOR * (1 - decay);
        }

        return { id: d.id, logOR, se, logSafe, seSafe: d.se_safety };
    });

    const maMort = runHKSJ(adjustedData, 'logOR', 'se');
    const maSafe = runHKSJ(adjustedData, 'logSafe', 'seSafe');

    const baseMort = 0.50;
    const baseBleed = 0.10;
    const trtMort = (baseMort * maMort.or) / (1 - baseMort + (baseMort * maMort.or));
    const trtBleed = (baseBleed * maSafe.or) / (1 - baseBleed + (baseBleed * maSafe.or));
    const livesSaved = Math.max(0, Math.round((baseMort - trtMort) * 100));
    const extraHarm = Math.max(0, Math.round((trtBleed - baseBleed) * 100));

    return { maMort, maSafe, livesSaved, extraHarm };
}

if (typeof module !== 'undefined' && module.exports) {
    module.exports = { runHKSJ, tInv975, runShockSim };
}
