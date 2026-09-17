#!/usr/bin/env python3
"""Exact checks for the RaQM gate-composition research note.

Python 3.10+, standard library only. No floating-point membership tests.
Run: python verify_gate_composition.py --output verification_results.json
Finite enumeration checks examples; the general results are proved in the note.
"""
from __future__ import annotations

import argparse
import json
from collections import deque
from fractions import Fraction as F
from itertools import permutations, product
from math import isqrt
from pathlib import Path


def conditional_probabilities(p):
    """Root probability and the conditional probabilities of occupied branches."""
    a = p[0] + p[1]
    return (a,) + ((p[0] / a,) if a else ()) + (
        (p[2] / (1 - a),) if a != 1 else ()
    )


def admissible(p, L, occupation=False):
    return (
        len(p) == 4 and all(x >= 0 for x in p) and sum(p) == 1
        and all((L * x).denominator == 1 for x in conditional_probabilities(p))
        and (not occupation or all((L * x).denominator == 1 for x in p))
    )


def probability_grid(L):
    result = set()
    for m, n, r in product(range(L + 1), repeat=3):
        a, b, c = F(m, L), F(n, L), F(r, L)
        result.add((a * b, a * (1 - b), (1 - a) * c, (1 - a) * (1 - c)))
    return result


def permute(v, pi):
    """pi[j] is the destination of input basis state j."""
    out = [0] * 4
    for j in range(4):
        out[pi[j]] = v[j]
    return tuple(out)


def tree_permutation(pi):
    return frozenset((frozenset(pi[:2]), frozenset(pi[2:]))) == frozenset(
        (frozenset((0, 1)), frozenset((2, 3)))
    )


def canonical(pi, d, L):
    """D_d P_pi normal form, with OUTPUT phase d[0] fixed to zero."""
    return tuple(pi), tuple((x - d[0]) % L for x in d)


def compose(left, right, L):
    """left after right, compared modulo one global phase."""
    pi, d = left
    sigma, e = right
    inv_pi = [pi.index(j) for j in range(4)]
    return canonical(
        [pi[sigma[j]] for j in range(4)],
        [d[j] + e[inv_pi[j]] for j in range(4)], L
    )


def inverse(g, L):
    pi, d = g
    return canonical([pi.index(j) for j in range(4)], [-d[pi[j]] for j in range(4)], L)


def basis_action(g):
    pi, d = g
    return [(pi[j], d[pi[j]]) for j in range(4)]


def check_composition_by_basis_action(left, right, L):
    l, r = basis_action(left), basis_action(right)
    sequential = [(l[r[j][0]][0], (r[j][1] + l[r[j][0]][1]) % L) for j in range(4)]
    direct = basis_action(compose(left, right, L))
    assert [x[0] for x in sequential] == [x[0] for x in direct]
    offsets = {(sequential[j][1] - direct[j][1]) % L for j in range(4)}
    assert len(offsets) == 1


def generated_group(L):
    ident = canonical(range(4), (0, 0, 0, 0), L)
    # X_A, CNOT(A->B), phase on A, phase on B, controlled phase.
    generators = [
        canonical((2, 3, 0, 1), (0, 0, 0, 0), L),
        canonical((0, 1, 3, 2), (0, 0, 0, 0), L),
        canonical(range(4), (0, 0, 1, 1), L),
        canonical(range(4), (0, 1, 0, 1), L),
        canonical(range(4), (0, 0, 0, 1), L),
    ]
    seen, todo = {ident}, deque([ident])
    while todo:
        g = todo.popleft()
        for h in generators:
            check_composition_by_basis_action(h, g, L)
            q = compose(h, g, L)
            if q not in seen:
                seen.add(q)
                todo.append(q)
    for g in seen:
        assert tree_permutation(g[0])
        assert compose(g, inverse(g, L), L) == ident
        assert compose(inverse(g, L), g, L) == ident
    assert len(seen) == 8 * L**3
    return len(seen)


def nearest_integer(x):
    """Exact nearest integer, ties away from zero; odd under sign reversal."""
    x = F(x)
    if x < 0:
        return -nearest_integer(-x)
    y = x + F(1, 2)
    return y.numerator // y.denominator


def primitive_diagonal_phases(a, b, c, L):
    """a,b,c are dimensionless energies in units epsilon_L, not radians.

    Return exponents of zeta for exp(-i H tau / hbar), quantising the
    physically specified local and interaction coefficients separately.
    """
    A, B, C = map(nearest_integer, (a, b, c))
    return (0, -B % L, -A % L, -(A + B + C) % L)


def entangling_circuit(L=8):
    # Exact amplitudes v[j] / sqrt(2)**h; all numerators are integers.
    v, h = (1, 0, 0, 0), 0
    rows = []

    def record(label):
        p = tuple(F(x*x, 2**h) for x in v)
        assert admissible(p, L, occupation=True)
        assert sum(p) == 1
        concurrence_squared = F(4 * (v[0]*v[3] - v[1]*v[2])**2, 2**(2*h))
        rows.append({"operation": label, "amplitude_numerators": list(v),
                     "sqrt2_denominator_power": h,
                     "probabilities": list(map(str, p)),
                     "concurrence_squared": str(concurrence_squared)})

    record("prepare |00>")
    for gate in ("H_A", "CNOT_A_to_B", "CZ", "CNOT_A_to_B", "H_A"):
        if gate == "H_A":
            v = (v[0]+v[2], v[1]+v[3], v[0]-v[2], v[1]-v[3])
            h += 1
        elif gate == "CNOT_A_to_B":
            v = permute(v, (0, 1, 3, 2))
        else:
            v = (v[0], v[1], v[2], -v[3])
        record(gate)
    assert rows[2]["concurrence_squared"] == "1"
    assert rows[-1]["probabilities"] == ["0", "0", "1", "0"]
    return rows


def verify():
    report = {"arithmetic": "exact integers and fractions; no numerical tolerances",
              "scope": "ray-level conditional grid and occupation-compatible subset; not a hidden-variable dynamics",
              "grid_checks": [], "monomial_group_checks": [], "counterexamples": {}}
    all_perms = list(permutations(range(4)))
    expected = {pi for pi in all_perms if tree_permutation(pi)}
    assert len(expected) == 8
    for L in (4, 8, 16):
        T = probability_grid(L)
        K = {p for p in T if admissible(p, L, occupation=True)}
        stabilisers = []
        for space in (T, K):
            good = {pi for pi in all_perms if all(permute(p, pi) in space for p in space)}
            assert good == expected
            stabilisers.append(len(good))
        assert len(T) == (L+1)*(L*L+1)
        report["grid_checks"].append({"L": L, "conditional_probability_points": len(T),
            "occupation_compatible_probability_points": len(K),
            "permutation_stabiliser_orders": stabilisers})
    for L in (4, 8):
        report["monomial_group_checks"].append({"L": L, "generated_order": generated_group(L),
            "expected_order": 8*L**3, "inverse_and_basis_action_checks": "passed"})

    witness = (F(1,4), F(1,4), F(1,2), F(0))
    swapped = permute(witness, (0,2,1,3))
    for L in (4,8,16,32,64):
        assert admissible(witness, L, occupation=True)
        assert not admissible(swapped, L)
    report["counterexamples"]["swap"] = {"input": list(map(str,witness)),
        "output": list(map(str,swapped)), "output_conditionals": list(map(str,conditional_probabilities(swapped))),
        "general_proof": "1/3 is not dyadic; valid for every L=2^M with M>=2"}
    p = (F(1,16), F(3,16), F(3,16), F(9,16))
    assert admissible(p,4) and not admissible(p,4,occupation=True)
    report["counterexamples"]["occupation"] = {"L":4,"probabilities":list(map(str,p))}

    # H_A maps the admissible (sqrt(3)/2)|00>+(1/2)|10> to a state
    # with p00=(2+sqrt(3))/4. Its polynomial has nonsquare discriminant.
    discriminant = (-16)**2-4*16
    assert isqrt(discriminant)**2 != discriminant
    report["counterexamples"]["hadamard_intermediate"] = {
        "probability": "(2+sqrt(3))/4", "irreducible_polynomial": "16*p^2-16*p+1",
        "discriminant": discriminant, "H_A_squared": "identity; endpoint membership does not ensure path membership"}

    L, a, b = 8, F(2,5), F(2,5)
    spectral = tuple(-nearest_integer(x) % L for x in (0,b,a,a+b))
    local = primitive_diagonal_phases(a,b,0,L)
    assert spectral == (0,0,0,L-1) and local == (0,0,0,0)
    report["counterexamples"]["independent_spectral_rounding"] = {
        "L":L,"local_coefficients_in_epsilon_L_units":[str(a),str(b)],
        "whole_spectrum_phase_exponents":list(spectral),"separately_quantised_exponents":list(local),
        "spurious_concurrence_on_plus_plus":"sin(pi/L)"}
    for a,b in product((F(-7,5),F(-1,2),F(2,5),F(3,2)), repeat=2):
        d=primitive_diagonal_phases(a,b,0,L)
        qa,qb=nearest_integer(a),nearest_integer(b)
        assert d==tuple(x%L for x in (0,-qb,-qa,-qa-qb))
    report["entangling_interference_circuit"] = entangling_circuit()
    report["status"] = "all checks passed"
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=Path("verification_results.json"))
    args=parser.parse_args()
    report=verify()
    args.output.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    for row in report["grid_checks"]:
        print(f"L={row['L']}: T probabilities={row['conditional_probability_points']}, "
              f"K probabilities={row['occupation_compatible_probability_points']}, "
              f"permutation orders={row['permutation_stabiliser_orders']}")
    for row in report["monomial_group_checks"]:
        print(f"L={row['L']}: generated projective monomial group order {row['generated_order']}")
    print(report["status"])
    print(f"Wrote {args.output}")


if __name__=="__main__":
    main()
