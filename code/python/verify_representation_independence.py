#!/usr/bin/env python3
"""Exact finite checks for RaQM representation independence, v0.1.
Python 3.10+, standard library only. Run this saved file to create a local
JSON record and assertion log beside it, under reproducibility/.
This verifies a finite Clifford/Bell model and rational identities, not the
all-resolution theorems by enumeration.
"""
from itertools import permutations, product
from fractions import Fraction as Q
from collections import Counter
from pathlib import Path
import datetime as dt
import hashlib
import json
import platform
import uuid


def transpose(a):
    return tuple(zip(*a))


def mul(a, b):
    return tuple(tuple(sum(x*y for x, y in zip(row, col))
                       for col in transpose(b)) for row in a)


def det(a):
    return (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
            -a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
            +a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))


def colors(values):
    seen = {}
    return tuple(seen.setdefault(v, len(seen)) for v in values)


def main():
    source = Path(__file__).resolve()
    now = lambda: dt.datetime.now(dt.timezone.utc).isoformat()
    start = now()
    record_dir = source.parent / 'reproducibility' / (
        'representation_python_' + dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
        + '_' + uuid.uuid4().hex[:8])
    record_dir.mkdir(parents=True)
    checks, summary = [], {}

    def check(test, label):
        checks.append({'label': label, 'passed': bool(test)})
        if not test:
            raise AssertionError(label)

    status, error = 'passed', None
    try:
        ident = ((1,0,0),(0,1,0),(0,0,1))
        reflect = ((1,0,0),(0,-1,0),(0,0,1))
        phase = ((0,-1,0),(1,0,0),(0,0,1))
        group = sorted(a for p in permutations(range(3)) for signs in product((-1,1), repeat=3)
                       if det(a := tuple(tuple(signs[i] if j == p[i] else 0 for j in range(3))
                                         for i in range(3))) == 1)
        index = {a:i for i,a in enumerate(group)}
        check(len(group) == 24, '24 distinct proper signed permutation rotations')
        tilde = lambda a: mul(mul(reflect, transpose(a)), reflect)
        for a in group:
            check(tilde(a) in index and tilde(tilde(a)) == a, 'unitary-transpose action is involutive')
            for b in group:
                check(mul(a,b) in index, 'exact group closure')
                check(tilde(mul(a,b)) == mul(tilde(b),tilde(a)), 'transpose reverses products')
        labels = tuple(product(range(24), repeat=2))
        position = {x:i for i,x in enumerate(labels)}
        obs = tuple(index[mul(group[a],tilde(group[b]))] for a,b in labels)
        check(len(set(obs)) == 24 and set(Counter(obs).values()) == {24}, '24 initial rays, each with 24 labels')
        ii, si = index[ident], index[phase]
        swap = tuple(si if i == ii else ii if i == si else i for i in range(24))
        fa = tuple(position[(swap[a], b)] for a,b in labels)
        fb = tuple(position[(a, swap[b])] for a,b in labels)
        for k in range(576):
            check(fa[fa[k]] == k and fb[fb[k]] == k, 'local updates are involutions')
            check(fa[fb[k]] == fb[fa[k]], 'local updates commute')
        x, y = position[(si,ii)], position[(ii,si)]
        check(obs[x] == obs[y] and obs[fa[x]] != obs[fa[y]], 'Bell aliases separate after one local swap')
        check(mul(transpose(group[obs[fa[x]]]),group[obs[fa[y]]]) == mul(phase,phase),
              'witness relative output rotation is phase squared')
        # Test all left multiplications, including non-dihedral exceptional gates.
        for v in group:
            f = tuple(index[mul(v,a)] for a in group)
            expected = tuple(index[mul(v,group[o])] for o in obs)
            actual = tuple(index[mul(group[f[a]],tilde(group[b]))] for a,b in labels)
            check(actual == expected, 'common left multiplication descends on all 576 Bell labels')
        initial = colors(obs)
        current, block_counts = initial, [len(set(initial))]
        while True:
            nxt = colors((current[k],current[fa[k]],current[fb[k]]) for k in range(576))
            if nxt == current:
                break
            current = nxt
            block_counts.append(len(set(current)))
        oracle = colors(tuple(obs[j] for j in (k,fa[k],fb[k],fa[fb[k]])) for k in range(576))
        check(current == oracle, 'partition refinement equals exhaustive four-word observation oracle')
        quotient_size = len(set(current))
        for f in (fa,fb):
            qmap = {}
            for k in range(576):
                qmap.setdefault(current[k],current[f[k]])
                check(qmap[current[k]] == current[f[k]], 'quotient transition is well defined')
            check(sorted(qmap.values()) == list(range(quotient_size)), 'quotient transition is a permutation')
        fibre_counts = [len({current[k] for k in range(576) if obs[k] == ray}) for ray in range(24)]
        for power in range(2,9):
            l = 2**power
            for dp in range(1,power+1):
                d = 2**dp
                z2 = 2 + l*sum((2*Q(n,d)-1)**2 for n in range(1,d))
                x2 = Q(l,2)*sum(1-(2*Q(n,d)-1)**2 for n in range(1,d))
                check(z2 == 2+Q(l*(d-1)*(d-2),3*d), 'exact axial second moment')
                check(x2 == Q(l*(d*d-1),3*d), 'exact transverse second moment')
                check(z2-x2 == 2-l+Q(l,d), 'exact second moment difference')
                check((z2 == x2) == (d == 2 and l == 4), 'unique sampled isotropic exception')
        summary.update(label_count=576, initial_ray_count=24,
                       refinement_block_counts=block_counts, quotient_class_count=quotient_size,
                       context_counts_by_ray=sorted(fibre_counts), maximum_context_count=max(fibre_counts),
                       conditional_fixed_width_bits=(max(fibre_counts)-1).bit_length(),
                       all_arithmetic='integer or rational; no floating-point checks')
    except Exception as exc:
        status, error = 'failed', repr(exc)
    record = dict(script=source.name, script_version='0.1', script_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                  started_utc=start, finished_utc=now(), python_version=platform.python_version(),
                  platform=platform.platform(), status=status, error=error,
                  passed_assertions=sum(c['passed'] for c in checks), total_assertions=len(checks),
                  summary=summary, checks=checks)
    (record_dir/'representation_independence_run_record.json').write_text(json.dumps(record,indent=2)+'\n')
    (record_dir/'representation_independence_assertion_log.txt').write_text(
        '\n'.join(('PASS: ' if c['passed'] else 'FAIL: ')+c['label'] for c in checks)+'\n')
    print(json.dumps({k:v for k,v in record.items() if k != 'checks'},indent=2))
    print('Local record:', record_dir)
    if status != 'passed':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
