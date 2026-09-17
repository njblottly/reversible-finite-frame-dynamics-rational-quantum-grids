#!/usr/bin/env python3
"""RaQM conditional context law, v0.1. Python 3.10+, standard library only.
Run this saved file; JSON and an assertion log are saved in reproducibility/.
Exact integer/rational checks in the 24-frame Clifford Bell submodel.
Uniform context weights and terminal Born responses are modelling assumptions.
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


def tr(a): return tuple(zip(*a))
def mm(a,b): return tuple(tuple(sum(x*y for x,y in zip(row,col)) for col in tr(b)) for row in a)
def det(a):
    return sum(a[0][j]*(a[1][(j+1)%3]*a[2][(j+2)%3]-a[1][(j+2)%3]*a[2][(j+1)%3]) for j in range(3))
def mean(ms): return tuple(tuple(sum(Q(m[i][j],len(ms)) for m in ms) for j in range(3)) for i in range(3))
def zero(): return ((0,0,0),(0,0,0),(0,0,0))
def overlap(r,s): return (1+sum(mm(tr(r),s)[i][i] for i in range(3)))/Q(4)


def main():
    path=Path(__file__).resolve()
    utc=lambda: dt.datetime.now(dt.timezone.utc).isoformat()
    started=utc(); checks=[]; summary={}; status='passed'; error=None
    folder=path.parent/'reproducibility'/('preparation_python_'+dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8])
    folder.mkdir(parents=True)
    def check(ok,label):
        checks.append(dict(label=label,passed=bool(ok)))
        if not ok: raise AssertionError(label)
    try:
        ident=((1,0,0),(0,1,0),(0,0,1)); jmat=((1,0,0),(0,-1,0),(0,0,1)); phase=((0,-1,0),(1,0,0),(0,0,1))
        group=sorted(a for p in permutations(range(3)) for s in product((-1,1),repeat=3)
                     if det(a:=tuple(tuple(s[i] if k==p[i] else 0 for k in range(3)) for i in range(3)))==1)
        check(len(set(group))==24,'24 distinct proper rotations')
        tilde=lambda a:mm(mm(jmat,tr(a)),jmat)
        swap=lambda a:phase if a==ident else ident if a==phase else a
        labels=tuple(product(group,repeat=2))
        obs=lambda x:mm(x[0],tilde(x[1]))
        action=lambda x,aa,bb:(swap(x[0]) if aa else x[0],swap(x[1]) if bb else x[1])
        fibres={r:tuple(x for x in labels if obs(x)==r) for r in group}
        signature=lambda x:tuple(obs(action(x,a,b)) for a,b in ((0,0),(1,0),(0,1),(1,1)))
        classes=Counter(signature(x) for x in labels)
        check(len(classes)==116,'116 retained context classes')
        for x in labels:
            check(signature(action(action(x,1,0),1,0))==signature(x),'exact context echo')
            for a,b in ((1,0),(0,1),(1,1)):
                check(classes[signature(x)]==classes[signature(action(x,a,b))],'counting measure on quotient is invariant')
        def output(r,aa=0,bb=0,safe=False):
            xs=fibres[r]
            if safe: xs=tuple(x for x in xs if all(a not in (ident,phase) for a in x))
            return mean([obs(action(x,aa,bb)) for x in xs])
        def prob(r,aa=0,bb=0,safe=False): return overlap(r,output(r,aa,bb,safe))
        for r in group:
            check(len(fibres[r])==24,'each pure preparation has 24 frame representatives')
            check(mean([obs(x) for x in fibres[r]])==r,'conditional law prepares its specified pure ray')
            one=output(r,1,0)
            expected=mean([r]*22+[mm(phase,r),mm(tr(phase),r)])
            check(one==expected,'one-sided uniform law equals the random-unitary channel')
            check(prob(r,1,0)==Q(23,24),'one-step Bell survival probability')
            safe_xs=[x for x in fibres[r] if all(a not in (ident,phase) for a in x)]
            check(len(safe_xs)>=20,'stationary-support preparation is nonempty')
            for a,b in ((0,0),(1,0),(0,1),(1,1)):
                check(output(r,a,b,True)==r,'stationary-support law reproduces identity for every command word')
            # Direct averaging over pushed quotient weights agrees with all 24 labels.
            local_classes=Counter(signature(x) for x in fibres[r])
            for pos,(a,b) in enumerate(((0,0),(1,0),(0,1),(1,1))):
                weighted=tuple(tuple(sum(Q(n,24)*key[pos][i][j] for key,n in local_classes.items()) for j in range(3)) for i in range(3))
                check(weighted==output(r,a,b),'induced quotient weights preserve terminal statistics')
        qphi=Counter(signature(x) for x in fibres[ident])
        check(sorted(qphi.values())==[1,1,1,21],'Bell Phi retained-context multiplicities')
        wrong=mean([key[1] for key in qphi])
        check(overlap(ident,wrong)==Q(3,4),'uniform quotient weights change the physical preparation law')
        # Fresh conditioning on each intermediate ray, rather than retained contexts.
        reset=mean([output(obs(action(x,1,0)),1,0) for x in fibres[ident]])
        check(overlap(ident,reset)==Q(265,288),'resampling destroys the exact two-step echo')
        pauli=[ident,((1,0,0),(0,-1,0),(0,0,-1)),((-1,0,0),(0,1,0),(0,0,-1)),((-1,0,0),(0,-1,0),(0,0,1))]
        ensembles=[pauli,[mm(phase,p) for p in pauli]]
        for en in ensembles:
            check(mean(en)==zero(),'both Bell-basis mixtures prepare I/4')
            check(mean([output(r,1,0) for r in en])==zero(),'one-sided ensemble consistency on the witness')
        t0=mean([output(r,1,1) for r in ensembles[0]])
        t1=mean([output(r,1,1) for r in ensembles[1]])
        expect=((0,Q(1,24),0),(Q(-1,24),0,0),(0,0,0))
        check(t0==expect and t1==tuple(tuple(-x for x in row) for row in expect),'exact two-sided mixture obstruction')
        # <X tensor Y> = -T_xy under the transpose convention for Bell correlations.
        p0=(1-t0[0][1])/2; p1=(1-t1[0][1])/2
        check((p0,p1)==(Q(23,48),Q(25,48)),'binary correlation probabilities separate equal-density preparations')
        # Finite fair-bit sampling: exact largest-remainder total variation formula.
        for m in range(1,33):
            for bits in range(0,9):
                n=2**bits; q,rem=divmod(n,m); counts=[q+1]*rem+[q]*(m-rem)
                tv=sum(abs(Q(c,n)-Q(1,m)) for c in counts)/2
                check(sum(counts)==n and tv==Q(rem*(m-rem),m*n),'exact balanced fair-bit approximation')
                check(tv<=Q(m,2*n),'finite-bit variation bound')
                check((tv==0)==(n%m==0),'exact uniformity iff fibre size divides bit-seed alphabet')
        summary=dict(label_count=576,retained_context_classes=116,
            phi_context_multiplicities=[1,1,1,21],one_step_survival='23/24',inverse_echo='1',
            two_steps_with_context_reset='265/288',uniform_quotient_one_step='3/4',
            equal_density_mixture_probabilities=['23/48','25/48'],
            mixture_probability_difference='1/24',arithmetic='integer and rational; no floating point',
            scope='Conditional predictions of specified preparation laws, not a unique RaQM prediction')
    except Exception as exc:
        status='failed'; error=repr(exc)
    record=dict(script=path.name,script_version='0.1',script_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                started_utc=started,finished_utc=utc(),python_version=platform.python_version(),platform=platform.platform(),
                status=status,error=error,passed_assertions=sum(c['passed'] for c in checks),total_assertions=len(checks),summary=summary,checks=checks)
    (folder/'preparation_law_run_record.json').write_text(json.dumps(record,indent=2)+'\n')
    (folder/'preparation_law_assertion_log.txt').write_text('\n'.join(('PASS: ' if c['passed'] else 'FAIL: ')+c['label'] for c in checks)+'\n')
    print(json.dumps({k:v for k,v in record.items() if k!='checks'},indent=2)); print('Local record:',folder)
    if status!='passed': raise SystemExit(1)

if __name__=='__main__': main()
