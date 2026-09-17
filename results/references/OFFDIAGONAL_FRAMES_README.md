# Off-diagonal marginals and basis/preparation selection

These files accompany version 0.5 of `main.tex`, dated 9 July 2026.

## What has been established

- An exact chord description of the admissible purifications of a fixed full-rank remote marginal, with an effective finite test when that marginal is specified algebraically.
- A restriction to polar or equatorial conditional rows when the marginal has nonzero coherence in the dyadic cyclotomic field.
- A fully classified family of off-diagonal marginals whose nonempty fixed-grid fibres have exactly `2 L` rays. A nominal local Hadamard has the sharp gap `Sqrt[2-Sqrt[2]]` for this family.
- An originally admissible state with irrational Schmidt probabilities, showing why diagonalising a marginal and imposing dyadic Schmidt weights adds a new condition.
- A bound on combined preparation, gate and active-basis tolerances while the remote marginal is kept exactly fixed.
- A proposed finite catalogue of local frames, with local updates, exact preservation of the remote marginal, commuting updates on distinct subsystems, and explicit approximation bounds for local-control histories.

The frame construction changes the admissible physical ray set. It is a candidate extension, not a proof that the original fixed grid already permits all local controls. The frame-label information budget, reversibility of the complete label update, consistency after forgetting labels, and a physical selection law involving Palmer's xi remain unresolved. The construction alone establishes neither preservation of the N_max prediction nor an experimental discrepancy law. Priority over the literature has not been established.

## Run the independent Python checks

Save `verify_offdiagonal_frames.py` in your local code folder. Use Python 3.10 or newer; no external packages or GPU are required.

```bash
python3 verify_offdiagonal_frames.py
```

The Python script also writes its own JSON and text record under a fresh directory:

```text
<script folder>/reproducibility/offdiagonal_python_<unique run identifier>/
    offdiagonal_frames_run_record.json
    offdiagonal_frames_assertion_log.txt
```

## Verification evidence delivered with this version

`offdiagonal_frames_authoring_record.json` records the actual Python execution in the authoring environment: **4,090 of 4,090 assertions passed**. 

The Python checks include exact enumeration at `L = 4, 8, 16`, using fractions, squarefree radicals and the cyclotomic power basis. They also include twelve 30-step numerical local-control histories, sampled preparation/gate/basis inequalities, and commuting-update checks. Numerical tolerance is `1e-10`; the largest sampled unitarity residual was below `1.1e-15`, and the largest remote-marginal residual below `8e-16`.

Exact enumeration verifies the stated finite instances. Floating-point residuals are numerical checks, not exact identities or interval certificates. The general mathematical arguments are the proofs in the manuscript.
