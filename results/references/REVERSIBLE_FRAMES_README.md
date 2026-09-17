# Reversible finite-frame control

This stage accompanies the single `main.tex`, version 0.6.1, dated 11 September 2026. 

## Mathematica repair: use script version 0.2

The supplied version 0.1 Mac run passed 901 assertions, then failed the first matched-permutation assertion. It completed no gate matching rows or local histories. The matcher incorrectly used `Return` inside `Do`, which exits that loop and lets the function reach its subsequent failure value. Version 0.2 replaces this control flow with explicit success flags. The mathematical tolerance has not been relaxed.

Before attempting gate matching, the repaired script checks identity, impossible and reassignment examples, and then all 512 bipartite graphs with three input and three output vertices against an exhaustive permutation oracle. It also logs the active catalogue and gate, and includes that context in any failure record.

The corrected matching core was executed in Mathics 10.0.1, an independent Wolfram Language interpreter. All 512 three-by-three graph cases agreed with exhaustive permutation search; identity, impossible and reassignment examples also returned the expected answers. The old fall-through behaviour was reproduced separately. This checks the repaired core, not the full numerical verifier in Mathematica.

## What the new results say

A permutation can replace the pointwise rounding update on the **same** finite frame catalogue. An explicit partition assigns equal Haar measure to every distinct frame, including the endpoint frames with collapsed phase coordinates. Hall's matching theorem then provides a permutation approximating each nominal left multiplication.

For dyadic `2 <= D` dividing `L`, the catalogue has

```text
N = 2 L + (D - 1) L^2
```

frames. A sufficient per-step projective operator error bound is

```text
B = 2 (Sqrt[2/(D - 1 + 2/L)] + 4 Sin[Pi/(2 L)]) .
```

Taking `D = L` retains an error of order `L^(-1/2)`. This conservative bound is uninformative on some small catalogues; numerical optimal matching costs can be much smaller. The proofs, rather than small numerical examples, establish its asymptotic behaviour.

The update is invertible on the frame labels, preserves the other subsystem's marginal, and commutes with updates on the other subsystem. It requires **no additional evolving state register**. The frame information already introduced in the previous stage is still required.

Inverse commands must use the inverse of the selected forward permutation. Independently matching the adjoint nominal matrix need not produce that inverse. General nominal gate products also need not be represented by products of their separately selected permutations. The manuscript explains why an exact nontrivial homomorphism from all of PU(2) into a finite permutation group is impossible.

If one instead insists on reproducing the original rounding exactly, a single clean-input implementation needs an auxiliary alphabet whose size is the largest preimage of that rounding map. This cost is necessary and sufficient for one step. It does not give an indefinitely reusable clean register.

The matching technique has classical precedents, credited in the manuscript. The explicit catalogue partition and resulting local-control bounds are the application developed here; priority over all literature has not been established. Label-independent physical dynamics, a selection law involving xi, operational ensemble consistency, and the physical frame/control information budget remain open. This result alone does not establish an unchanged N_max prediction.

## Evidence supplied with this version

`reversible_frames_authoring_record.json` records the actual Python run here: **1,612 of 1,612 assertions passed**. It includes all 288 finite maps on sets of sizes one through four, exact rational partition checks, twelve numerical gate matchings, and four 24-step histories followed by exact inverse label updates. Physical return residuals were below `2e-15`.
