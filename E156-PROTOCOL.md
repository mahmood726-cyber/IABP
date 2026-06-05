# E156-PROTOCOL — ShockLogic (Mechanical Support in Cardiogenic Shock)

- **Project:** IABP (GitHub repo `IABP`, user `mahmood726-cyber`)
- **Revived:** 2026-06-05 (from a single-file `IABP.html` dump)
- **Type:** single-file offline browser tool + Node-testable engine
- **Dashboard:** GitHub Pages (`index.html`)

## What changed in the revival

- Made **fully offline** (removed the Google Fonts CDN link).
- Extracted the HKSJ core into a pure `engine.js` (single source of truth).
- **Fixed two statistical correctness issues:** added the missing HKSJ variance
  floor (was anti-conservative when Q<k−1) and replaced a hardcoded t-value with
  `qt(0.975, k−1)`.
- Removed the Blob web-worker; the page now calls `engine.js` synchronously.
- Added `tests.js` (21 assertions) and the Pages scaffold; renamed to
  `index.html`.

## Body (E156 draft — CURRENT BODY)

Which mechanical circulatory support strategy improves survival in cardiogenic
shock without trading mortality for catastrophic bleeding? This dashboard pools
mortality and safety odds ratios for Impella, VA-ECMO and IABP from the named
shock trials and projects the net effect onto a 100-patient cohort as lives
saved versus extra harm. Pooling uses a Hartung–Knapp–Sidik–Jonkman
random-effects model with DerSimonian–Laird τ², and the user can stress the
result by STEMI subgroup and by time-to-support. A revival audit corrected two
anti-conservative defects — a missing HKSJ variance floor and a hardcoded
critical value in place of `qt(0.975, k−1)` — so the intervals now widen
honestly with few studies. Under the corrected model no device shows a robust
mortality benefit, and the active devices carry a clear excess of major
bleeding and vascular harm, so the net competing-risk projection is close to
neutral. The honest read is that current shock-support evidence does not
support routine device escalation for survival, and the wide intervals reflect
genuine uncertainty rather than analytic precision. The tool is a transparent
sensitivity aid, not a treatment rule.

SUBMITTED: [ ]
