# ShockLogic — Mechanical Support in Cardiogenic Shock

A single-file, **fully offline** dashboard that pools device odds ratios
(Impella/DanGer, VA-ECMO/ECLS-SHOCK, IABP-SHOCK II, ISAR-SHOCK, IMPRESS) for
mortality and major-bleeding/vascular harm, with subgroup (STEMI) and
time-to-support sensitivity, a forest plot, and a competing-risk waffle
("lives saved vs extra harm" per 100 patients).

**Live app:** open `index.html` (or the GitHub Pages link). No build, no network,
no external CDN.

## Layout

```
index.html   single-file UI (loads engine.js, calls it synchronously)
engine.js    pure HKSJ meta-analysis core — runs in Node and the browser
tests.js     Node test harness, 21 assertions
LICENSE      MIT
```

## Statistical core (`engine.js`)

`runHKSJ(studies, effKey, seKey)` pools log odds ratios with a
**Hartung–Knapp–Sidik–Jonkman** model: DerSimonian–Laird τ², random-effects
weighted mean, and the HK variance for the confidence interval.
`runShockSim({devices, subgroup, timeDelay})` runs the full pipeline (context
adjustment → pooling → baseline-risk competing-risk projection).

## Fixes applied during revival (2026-06-05)

1. **HKSJ variance floor.** The original divided the HK statistic with no floor,
   so when `Q < k−1` the interval narrowed *below* the random-effects interval
   (anti-conservative). Now floored at the classical RE variance
   (`max(1, q)·Var_RE`), per Knapp–Hartung / IntHout.
2. **Proper t critical value.** The original hardcoded 2.78 (k≤5) / 2.26 (else).
   Now uses `qt(0.975, k−1)` via an exact table (df 1–30) + Cornish–Fisher
   beyond — e.g. k=2 → 12.71 (was 2.78), k=3 → 4.30.
3. **Worker removed.** The Blob web-worker (a `file://` fragility for a trivial
   computation) was replaced by a synchronous call into `engine.js`.
4. **Offline.** Google Fonts CDN link removed; the app loads no external
   resource.

## Tests

```
node tests.js
# 21 passed, 0 failed
```

Includes t-quantiles vs R `qt()`, a floor-engaging identical-studies case
(seHKSJ must equal `√(1/200)`, not collapse to 0), and a hand-computed
heterogeneous case (τ²≈0.367, pooled OR≈0.769, I²≈95.8%).

## Caveats

HKSJ is preferred over DL-Wald for small *k* but with very few studies the
`t_{k−1}` intervals are wide by design (k=2 ⇒ t≈12.7). The bundled effect sizes
are illustrative estimates drawn from the named shock trials, not a registered
extraction. Hypothesis-generating only. MIT licensed.
