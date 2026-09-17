#!/usr/bin/env python3
"""Exact data and editable LaTeX for the RaQM paper's probability-history figure.

Python 3.10+, standard library only. Run in PyCharm or:
    python3 reproduce_paper_figures.py
Outputs are written under paper_figures/ beside this saved script.
Optional: --output-dir PATH. Add --latex-fragment only if you also want a
plot-only LaTeX fragment. No large matching search or old verifier is needed.
The self-contained manuscript already embeds the figure; running this script
is optional and provides a local record of its exact rational source data.
"""
from fractions import Fraction
from pathlib import Path
from datetime import datetime, timezone
import argparse
import csv
import hashlib
import json


class F:
    """Exact Q(sqrt(2),i); basis masks 0:1, 1:sqrt(2), 2:i, 3:i sqrt(2)."""
    def __init__(self, x=0):
        self.c = ({k: Fraction(v) for k, v in x.items() if v}
                  if isinstance(x, dict) else ({0: Fraction(x)} if x else {}))

    def __add__(self, other):
        other = field(other)
        out = self.c.copy()
        for k, v in other.c.items():
            out[k] = out.get(k, Fraction(0)) + v
        return F(out)

    __radd__ = __add__

    def __neg__(self):
        return F({k: -v for k, v in self.c.items()})

    def __sub__(self, other):
        return self + -field(other)

    def __mul__(self, other):
        other = field(other)
        out = {}
        for a, x in self.c.items():
            for b, y in other.c.items():
                common = a & b
                factor = (2 if common & 1 else 1) * (-1 if common & 2 else 1)
                key = a ^ b
                out[key] = out.get(key, Fraction(0)) + factor * x * y
        return F(out)

    __rmul__ = __mul__

    def conj(self):
        return F({k: -v if k & 2 else v for k, v in self.c.items()})

    def rational(self):
        if set(self.c) - {0}:
            raise ValueError('Expected a rational probability, got ' + str(self.c))
        return self.c.get(0, Fraction(0))


def field(x):
    return x if isinstance(x, F) else F(x)


def matrix(rows):
    return tuple(tuple(field(x) for x in row) for row in rows)


def mm(a, b):
    return tuple(tuple(sum((x * y for x, y in zip(row, col)), F())
                       for col in zip(*b)) for row in a)


def mv(a, v):
    return tuple(sum((x * y for x, y in zip(row, v)), F()) for row in a)


def kron(a, b):
    return tuple(tuple(a[i // len(b)][j // len(b[0])] * b[i % len(b)][j % len(b[0])]
                       for j in range(len(a[0]) * len(b[0])))
                 for i in range(len(a) * len(b)))


def probability_a1(v):
    return sum((v[j].conj() * v[j] for j in (2, 3)), F()).rational()


def calculate(horizon=16):
    halfroot = F({1: Fraction(1, 2)})
    imag = F({2: 1})
    I = matrix([[1, 0], [0, 1]])
    H = matrix([[halfroot, halfroot], [halfroot, -halfroot]])
    T = matrix([[1, 0], [0, halfroot * (1 + imag)]])
    S = mm(T, T)
    Z = matrix([[1, 0], [0, -1]])
    C = matrix([[1, 0, 0, 0], [0, 1, 0, 0],
                [0, 0, 0, 1], [0, 0, 1, 0]])
    HA, TA, SA = kron(H, I), kron(T, I), kron(S, I)
    RB = kron(I, mm(Z, H))
    ideal_step = mm(C, mm(HA, TA))
    original_odd = mm(C, HA)
    rematched_odd = mm(C, mm(RB, HA))
    even = mm(C, mm(HA, SA))
    ideal = original = rematched = tuple(map(F, (1, 0, 0, 0)))
    rows = []
    for n in range(horizon + 1):
        # Controller vectors are inner Clifford-frame vectors G_n|00>.
        # The omitted outer T_A^(n mod 2) is diagonal and cannot change A=1.
        p, q, r = map(probability_a1, (ideal, original, rematched))
        rows.append(dict(n=n, ideal=p, original=q, rematched=r,
                         original_gap=abs(q-p), rematched_gap=abs(r-p)))
        odd_next = (n + 1) % 2
        ideal = mv(ideal_step, ideal)
        original = mv(original_odd if odd_next else even, original)
        rematched = mv(rematched_odd if odd_next else even, rematched)
    return rows


def history_figure(rows):
    def coords(key):
        return ' '.join(f"({r['n']},{float(r[key]):.16g})" for r in rows)
    return r'''% Generated from exact rational probabilities by reproduce_paper_figures.py.
\begin{figure}[tbp]
\centering
\begin{tikzpicture}
\begin{groupplot}[
 group style={group size=1 by 2,vertical sep=1.15cm},
 width=0.92\linewidth,height=5.0cm,
 xmin=-0.3,xmax=16.3,xtick={0,2,4,6,8,10,12,14,16},
 tick label style={font=\small},label style={font=\small},
 title style={font=\small},grid=major,
 grid style={gray!16},axis line style={gray!65},
 legend style={font=\footnotesize,draw=none,fill=white},
 clip=true]
\nextgroupplot[ymin=-0.04,ymax=1.08,ylabel={$P(A=1)$},
 title={(a) Exact terminal probabilities},
 legend style={at={(0.5,1.23)},anchor=south},legend columns=3]
\addplot[black,thick,mark=*,mark size=1.6pt] coordinates {IDEAL};
\addlegendentry{Ideal commands}
\addplot[raqmOrange,thick,dashed,mark=square*,mark size=2pt]
 coordinates {ORIGINAL};
\addlegendentry{Original CNOT table}
\addplot[raqmBlue,thick,dashdotted,mark=o,mark size=2.6pt]
 coordinates {REMATCHED};
\addlegendentry{Rematched CNOT table}
\nextgroupplot[ymin=-0.025,ymax=0.70,ylabel={Absolute probability error},
 xlabel={Repeated blocks $n$ (three elementary commands per block)},
 title={(b) Discrepancy from the ideal history}]
\addplot[gray,densely dashed,forget plot] coordinates {(7,0) (7,0.67)};
\addplot[gray,densely dashed,forget plot] coordinates {(15,0) (15,0.67)};
\addplot[black,dotted,thick,forget plot] coordinates {(0,0.3333333333333333) (16,0.3333333333333333)};
\node[font=\footnotesize,anchor=south west,fill=white,inner sep=1pt]
 at (axis cs:0.3,0.3333333333333333) {$\delta=1/3$};
\addplot[raqmOrange,thick,dashed,mark=square*,mark size=2pt]
 coordinates {OGAP};
\addplot[raqmBlue,thick,dashdotted,mark=o,mark size=2.6pt]
 coordinates {RGAP};
\node[font=\footnotesize,anchor=south,fill=white,inner sep=1pt]
 at (axis cs:7,0.64) {$n=7$};
\node[font=\footnotesize,anchor=south,fill=white,inner sep=1pt]
 at (axis cs:15,0.64) {$n=15$};
\end{groupplot}
\end{tikzpicture}
\caption{Exact histories of the repeated block $(T_A,H_A,C)$ from
$\ket{00}$, for the specified original table with $K_0=K_1=H$ and the
CNOT rematching of \eqref{eq:cnot-countertable}. Both use the same label
count, calibration anchors, Hadamard rule and coarse error bounds.
The first strict $1/3$ crossing is at seven blocks for this original
table and fifteen for the rematched table. At seven, their errors are
$7/16$ and $1/16$, respectively. Every point holds for every fixed prior
on the specified compensated phase orbit. Points represent integer
depths; connecting lines guide the eye. These are two particular
controllers, not an optimisation over matchings or an asymptotic
refinement result.}
\label{fig:cnot-histories}
\end{figure}
'''.replace('IDEAL', coords('ideal')).replace('ORIGINAL', coords('original')).replace(
        'REMATCHED', coords('rematched')).replace('OGAP', coords('original_gap')).replace(
        'RGAP', coords('rematched_gap'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path,
                        default=Path(__file__).resolve().parent / 'paper_figures')
    parser.add_argument('--latex-fragment', action='store_true',
                        help='Also write a plot-only LaTeX fragment, not a manuscript.')
    args = parser.parse_args()
    rows = calculate()
    threshold = Fraction(1, 3)
    first = lambda key: next(r['n'] for r in rows if r[key] > threshold)
    checks = {
        'original_first_strict_crossing_is_7': first('original_gap') == 7,
        'rematched_first_strict_crossing_is_15': first('rematched_gap') == 15,
        'seventh_block_ideal_is_1_over_16': rows[7]['ideal'] == Fraction(1, 16),
        'seventh_block_original_gap_is_7_over_16': rows[7]['original_gap'] == Fraction(7, 16),
        'seventh_block_rematched_gap_is_1_over_16': rows[7]['rematched_gap'] == Fraction(1, 16),
        'fifteenth_block_rematched_gap_is_95_over_256': rows[15]['rematched_gap'] == Fraction(95, 256),
        'all_probabilities_are_in_unit_interval': all(0 <= r[k] <= 1 for r in rows
                                                        for k in ('ideal','original','rematched')),
    }
    if not all(checks.values()):
        raise AssertionError(checks)
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    serial = [{k: v if k == 'n' else str(v) for k, v in row.items()} for row in rows]
    with (out / 'raqm_history_data.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(serial)
    record = {
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'status': 'passed', 'arithmetic': 'exact Q(sqrt(2),i), rational probabilities',
        'scope': 'Two specified controllers; not a full matching verification or optimisation.',
        'checks': checks, 'rows': serial,
    }
    (out / 'paper_figure_run_record.json').write_text(json.dumps(record, indent=2) + '\n')
    if args.latex_fragment:
        (out / 'raqm_history_figure.tex').write_text(history_figure(rows))
    print('Figure data checks: passed')
    print('Output:', out)
    print('The manuscript already contains the figure. No Overleaf image upload is required.')


if __name__ == '__main__':
    main()
