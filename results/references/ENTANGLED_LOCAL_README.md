# Entangled local-control verification

The manuscript is still the single file `main.tex` (working version
0.4). The earlier verification scripts remain separate records of
earlier results; these new checks concern the entangled local-control section.

## What was established

- Complete local-fibre classifications for diagonal unchanged marginals,
  under both of the manuscript's specified grids.
- A sharp local approximation gap of `2 sin(pi/16)` for nonmaximal,
  computational-basis Schmidt seeds when acting on the root qubit. It persists
  for arbitrarily large dyadic grids, including the conditional grid.
- A constructive rule for histories on the other qubit with a fixed diagonal
  marginal, including nonmaximal entanglement.
- A constructive rule for alternating local controls on either qubit throughout
  the admissible maximally entangled sector, with a proved per-step upper bound.

The classification is chart-specific. General off-diagonal marginals,
preparation/chart selection, and the hidden-variable outcome law remain open.
The existing fixed-marginal existence criterion is standard purification
uniqueness, now explicitly credited to Watrous, Theorem 2.12. Literature priority
for the grid-specific results has not been established by the preliminary search.

## Current verification status

- Python: passed in the authoring environment, 10,546/10,546 assertions.
  `entangled_local_authoring_record.json` is that run.
- Exact Python checks: 100 full probability-grid fibre classifications at
  L=4,8,16; phase enumeration at L=4; 3,844 dyadic valuation cases.
- Numerical Python checks: unitary approximation, sharp gap, 40-step alternating
  Bell histories and 30-step B-only histories, tolerance 1e-10, seed 20260908.
- Mathematica: record confirms 1,859/1,859 assertions
  passed in Mathematica 13.3.0 on MacOSX-ARM64, using 50-digit precision and
  tolerance 1e-25. All 24 fibre rows match the Python results; all 15 coverage
  rows satisfy their bounds. Eight 20-step local histories also passed.
  Execution was interactive, so the executed source hash was not captured.
- Numerical checks are not interval certificates. General proofs are in the
  manuscript. Neither finite enumeration nor a zero reported numerical residual
  establishes the result for all L on its own.
