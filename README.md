# Reversible Finite-Frame Dynamics for Rational Quantum Grids

Exact checks for the rational state grids and gate constructions studied in the paper. The gate-composition tests distinguish the conditional-amplitude grid $T_L$ from the occupation-compatible subset $K_L$, check their monomial symmetries, and evaluate the finite circuit examples.

[Paper — DOI: 10.5281/zenodo.22313144](https://doi.org/10.5281/zenodo.22313144)

## Manuscript


The `main.tex` and `references.bib` were produced in an Overleaf project. The manuscript uses APS REVTeX 4.2; for compatibility, both figures are embedded in the LaTeX source.

## Gate-composition checks

Python 3.10 or newer, standard library only. Run from the directory containing `verify_gate_composition.py`:

```bash
python3 verify_gate_composition.py --output verification_results.json
```

In PyCharm, select a Python 3.10+ interpreter and run the same script. Without arguments, it writes `verification_results.json` in the working directory. An existing file at that path is overwritten.

The script uses exact fractions for grid membership. It tests all 24 basis permutations at $L=4,8,16$, enumerates the 512- and 4096-element projective monomial groups at $L=4,8$, and checks composition, inverses, circuit probabilities and the finite counterexamples.

Expected counts:

| $L$ | $T_L$ probability vectors | $K_L$ probability vectors | Permutations preserving $T_L$ | Permutations preserving $K_L$ |
|---:|---:|---:|---:|---:|
| 4 | 85 | 27 | 8 | 8 |
| 8 | 585 | 77 | 8 | 8 |
| 16 | 4369 | 233 | 8 | 8 |

The permutation counts exclude diagonal phase choices. The full projective monomial group has order $8L^3$.

For the Mathematica checks, no add-on packages are required. Save a notebook beside `verify_gate_composition.wl`, then evaluate:

```wolfram
SetDirectory[NotebookDirectory[]];
Get["verify_gate_composition.wl"]
```

With `wolframscript` installed, the equivalent command is:

```bash
wolframscript -file verify_gate_composition.wl
```

The script prints the results and aborts on a failed assertion. These finite checks support the proofs; they do not enumerate the state spaces at large $L$.

The commands above cover gate composition, not the complete verification suite.

## Figure data

Run from the directory containing the figure script:

```bash
python3 reproduce_paper_figures.py
```

It writes `paper_figures/raqm_history_data.csv` and `paper_figures/paper_figure_run_record.json`. The CSV supplies the exact history data for the circuit comparison. The schematic matching diagram is defined directly in the manuscript.

## Scope

The calculations concern explicitly defined ray grids and finite controllers. They do not establish a complete encoding of Palmer’s hidden variable $\xi$, unrestricted matching optimality, or new experimental bounds.

Palmer’s source papers are [arXiv:2510.02877](https://arxiv.org/abs/2510.02877) and [arXiv:2602.16382](https://arxiv.org/abs/2602.16382). Their PDFs are not required to run the checks.

## Citation

```bibtex
@misc{Bruzzese2026FiniteFrame,
  author = {Bruzzese, Josef},
  title  = {Reversible Finite-Frame Dynamics for Rational Quantum Grids:
            Preparation Laws and Matching-Dependent Statistical Predictions},
  year   = {2026},
  doi    = {10.5281/zenodo.22313144},
  url    = {https://doi.org/10.5281/zenodo.22313144}
}
```