/*
 * Node tests for the ShockLogic engine. Run: node tests.js
 * Expectations hand-computed / from R qt() — independent of the engine.
 */
const { runHKSJ, tInv975, runShockSim } = require('./engine.js');

let pass = 0, fail = 0;
function ok(name, cond, detail) {
    if (cond) { pass++; console.log('  ok  ' + name); }
    else { fail++; console.log(' FAIL ' + name + (detail ? '  -> ' + detail : '')); }
}
const close = (a, b, tol) => Math.abs(a - b) < (tol || 1e-3);

// --- tInv975 vs R qt(0.975, df) ---
ok('t df=1  == 12.706', close(tInv975(1), 12.706));
ok('t df=2  == 4.303', close(tInv975(2), 4.303));
ok('t df=4  == 2.776', close(tInv975(4), 2.776));
ok('t df=9  == 2.262', close(tInv975(9), 2.262));
ok('t df=30 == 2.042', close(tInv975(30), 2.042));
ok('t df=100 ~ 1.984 (Cornish-Fisher)', close(tInv975(100), 1.98397, 2e-3), 'got ' + tInv975(100));
ok('t monotone decreasing', tInv975(3) > tInv975(10) && tInv975(10) > tInv975(60));

// --- single study passthrough ---
const s1 = runHKSJ([{ id: 'A', logOR: Math.log(0.8), se: 0.1 }], 'logOR', 'se');
ok('k=1: or == 0.8', close(s1.or, 0.8));
ok('k=1: I2 == 0', s1.I2 === 0);

// --- HKSJ FLOOR: two identical studies => qStat=0, floored to RE variance ---
// muRE=ln(0.8), wSum=200, varRE=1/200, seHKSJ=sqrt(0.005)=0.070711, t=12.706
const ident = runHKSJ([
    { id: 'A', logOR: Math.log(0.8), se: 0.1 },
    { id: 'B', logOR: Math.log(0.8), se: 0.1 }
], 'logOR', 'se');
ok('floor: or == 0.8', close(ident.or, 0.8));
ok('floor: tau2 == 0', close(ident.tau2, 0, 1e-9));
const seFloored = Math.log(ident.hi / ident.or) / tInv975(1); // recover seHKSJ from hi
ok('floor: seHKSJ == sqrt(1/200) (NOT 0)', close(seFloored, Math.sqrt(1 / 200), 1e-4), 'got ' + seFloored);
ok('floor: CI has width (not degenerate)', ident.hi > ident.or && ident.lo < ident.or);

// --- heterogeneous two-study case: hand-computed ---
// logOR1=ln(0.5) se=0.1; logOR2=ln(1.2) se=0.15
// tau2=0.36695, or(muRE)=0.76910, I2=95.76%, qStat=1.00009, seHKSJ=0.43768
const het = runHKSJ([
    { id: 'A', logOR: Math.log(0.5), se: 0.10 },
    { id: 'B', logOR: Math.log(1.2), se: 0.15 }
], 'logOR', 'se');
ok('het: tau2 ~ 0.36695', close(het.tau2, 0.36695, 2e-3), 'got ' + het.tau2);
ok('het: or ~ 0.76910', close(het.or, 0.76910, 2e-3), 'got ' + het.or);
ok('het: I2 ~ 95.76%', close(het.I2, 95.76, 0.1), 'got ' + het.I2);

// --- runShockSim integration ---
const DB = [
    { id: "DanGer (Impella)", or: 0.74, se: 0.12, or_safety: 2.50, se_safety: 0.20, active: true },
    { id: "ECLS-SHOCK (ECMO)", or: 0.98, se: 0.10, or_safety: 1.80, se_safety: 0.15, active: true },
    { id: "IABP-SHOCK II", or: 0.96, se: 0.08, or_safety: 1.05, se_safety: 0.10, active: true }
];
const sim = runShockSim({ devices: DB, subgroup: 'all', timeDelay: 0 });
ok('sim: returns results object', sim && sim.maMort && sim.maSafe);
ok('sim: mortality OR in (0,2)', sim.maMort.or > 0 && sim.maMort.or < 2);
ok('sim: livesSaved is a non-negative integer', Number.isInteger(sim.livesSaved) && sim.livesSaved >= 0);
ok('sim: safety OR > 1 (devices add harm)', sim.maSafe.or > 1);
ok('sim: no active devices -> null', runShockSim({ devices: DB.map(d => ({ ...d, active: false })), subgroup: 'all', timeDelay: 0 }) === null);

console.log('\n' + pass + ' passed, ' + fail + ' failed');
process.exit(fail === 0 ? 0 : 1);
