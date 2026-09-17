# RaQM: consolidated manuscript and local reproducibility

## Verification status

- The original `verify_gate_composition.wl` was successfully run again on
  7 September 2026. 
  - `raqm_gate_composition.tex` version 0.2 now records that success. Its
  mathematics is unchanged from version 0.1. No PDF was regenerated.
- The original Python script and the new `verify_nominal_gates.py` passed in
  the authoring environment. A run on your Mac is still required to create
  your Mac-local record.
- The 16 September 2026 `verify_nominal_gates.wl` reports
  that all mathematical checks passed, but automatic local-record export
  failed. The updated historical JSON records both facts.
- Logging revision 0.2 fixes filename resolution and the `AssociateTo` call,
  and checks whether both exports succeeded. Its mathematical assertion block
  is unchanged. Your subsequent exported JSON and summary log confirm that
  revision 0.2 saved the records successfully: all 1,131 assertions passed
  in Mathematica 13.3.0 on MacOSX-ARM64. The run used interactive evaluation,
  so no executed source hash was captured. No rerun is needed to establish
  the successful checks and export already recorded.

## The consolidated manuscript

1. An exact precision-allocation decomposition of the occupation-compatible grid.
2. An explicit quantiser with uniform sufficient covering bounds.
3. A smallest-operator-distance unitary correction for a selected output ray.
4. An accumulated error bound for input-dependent selected gate sequences.
5. A counterexample to a uniformly small strictly local repair on the full
   occupation-compatible grid.
6. A constructive local-control rule for product states with reserved precision.
7. A fixed-marginal criterion for the remaining entangled local-control problem.

The unrestricted correction may act on both qubits even when the nominal
gate acts on one. The local product-sector construction avoids that issue
under a stated preparation restriction. Neither result yet supplies Palmer's
hidden-variable setting-selection law or a state-independent quantum channel.

The new Python program separates exact fraction checks from floating-point
matrix checks. Its matrix tolerance is 1e-10, and its random sample seed is
20260907. The Mathematica companion uses 50-digit working precision and
tolerance 1e-25. Numerical samples check the construction; the general proofs
are in the consolidated manuscript, and no rate is claimed to be optimal.

Wolfram documentation for the record format and source hashes:

- [FileHash](https://reference.wolfram.com/language/ref/FileHash.html)
- [RawJSON](https://reference.wolfram.com/language/ref/format/RawJSON.html)
