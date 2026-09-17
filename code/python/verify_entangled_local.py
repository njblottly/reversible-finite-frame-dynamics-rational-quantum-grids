#!/usr/bin/env python3
"""Standalone standard-library verification of the new entangled local-control results.
Run: python3 verify_entangled_local.py
Writes fresh local JSON and text records beside this script. Exact enumeration
and numerical checks are separate; numerical checks are not interval proofs.
"""
from fractions import Fraction as F
from pathlib import Path
from datetime import datetime, timezone
from itertools import product
import cmath
import hashlib
import json
import math
import platform
import random
import sys
import traceback
import uuid

TOL = 1e-10
SEED = 20260908
checks = []

def check(test, label):
    checks.append({'label': label, 'passed': bool(test)})
    if not test:
        raise AssertionError(label)

def grid(L, occupation=False):
    out = set()
    for m, n, r in product(range(L+1), repeat=3):
        a,b,c = F(m,L),F(n,L),F(r,L)
        p=(a*b,a*(1-b),(1-a)*c,(1-a)*(1-c))
        if not occupation or all((L*x).denominator==1 for x in p):
            out.add(p)
    return out

def fibre_probabilities(g,q,side):
    if side=='A':
        return [p for p in g if p[0]+p[2]==q and p[0]*p[1]==p[2]*p[3]]
    return [p for p in g if p[0]+p[1]==q and p[0]*p[2]==p[1]*p[3]]

def expected_t(L,q,side,kind):
    if side=='A' and q != F(1,2):
        return {F(0),F(1)} | ({F(1,2)} if kind=='T' or L%(2*q.denominator)==0 else set())
    D=L if kind=='T' else (L//2 if side=='A' else L//q.denominator)
    return {F(i,D) for i in range(D+1)}

def phase_count_exhaustive(ps,L,side):
    """Count all relative phases by exact modular cancellation of equal magnitudes."""
    count=0
    pairs=((0,1),(2,3)) if side=='A' else ((0,2),(1,3))
    for p in ps:
        occupied=[j for j in range(4) if p[j]]
        for choices in product(range(L),repeat=len(occupied)-1):
            phase=dict(zip(occupied,[0,*choices]))
            (i,j),(k,l)=pairs
            if p[i]*p[j]==0:
                ok=p[k]*p[l]==0
            else:
                ok=(phase[j]-phase[i]-phase[l]+phase[k])%L==L//2
            count+=bool(ok)
    return count

def mm(a,b):
    return [[sum(a[i][k]*b[k][j] for k in range(2)) for j in range(2)] for i in range(2)]

def adj(a):return [[a[j][i].conjugate() for j in range(2)] for i in range(2)]
def trans(a):return [list(x) for x in zip(*a)]
def scale(a,z):return [[z*x for x in row] for row in a]
def sub(a,b):return [[a[i][j]-b[i][j] for j in range(2)] for i in range(2)]
def frob(a):return math.sqrt(sum(abs(x)**2 for row in a for x in row))
def opnorm(a):
    g=mm(adj(a),a)
    return math.sqrt(max(0,(g[0][0].real+g[1][1].real+math.hypot(g[0][0].real-g[1][1].real,2*abs(g[0][1])))/2))
I=[[1+0j,0j],[0j,1+0j]]

def param(t,b,g):
    x,y=math.sqrt(t),math.sqrt(1-t)
    return [[x+0j,y*cmath.exp(1j*b)],[y*cmath.exp(1j*g),-x*cmath.exp(1j*(b+g))]]

def quantise(u,D,L):
    t=min(1,max(0,abs(u[0][0])**2))
    if abs(u[0][0]) < 1e-14:
        t=0.; lam=cmath.phase(u[0][1]); beta=0.; gamma=cmath.phase(u[1][0])-lam
    elif abs(u[0][1]) < 1e-14:
        t=1.; lam=cmath.phase(u[0][0]); beta=0.; gamma=cmath.phase(-u[1][1])-lam
    else:
        lam=cmath.phase(u[0][0]); beta=cmath.phase(u[0][1])-lam; gamma=cmath.phase(u[1][0])-lam
    n=max(0,min(D,round(D*t))); b=round(beta*L/(2*math.pi))%L; g=round(gamma*L/(2*math.pi))%L
    cert={'t':F(n,D),'beta_index':b,'gamma_index':g}
    return scale(param(float(cert['t']),2*math.pi*b/L,2*math.pi*g/L),cmath.exp(1j*lam)),cert

def random_unitary(rng):
    return scale(param(rng.random(),rng.uniform(-math.pi,math.pi),rng.uniform(-math.pi,math.pi)),cmath.exp(1j*rng.uniform(-math.pi,math.pi)))

def seed_matrix(q):return [[math.sqrt(float(q))+0j,0j],[0j,math.sqrt(float(1-q))+0j]]
def chord(a,b):
    z=sum(x.conjugate()*y for ar,br in zip(a,b) for x,y in zip(ar,br))
    aligned=scale(b,z.conjugate()/abs(z)) if abs(z)>0 else b
    return frob(sub(a,aligned))

def run():
    rows=[]; classifications=0
    for L in (4,8,16):
        for kind in ('T','K'):
            g=grid(L,kind=='K')
            for q in (F(n,L) for n in range(1,L)):
                for side in ('A','B'):
                    ps=fibre_probabilities(g,q,side)
                    ts={p[0]/q for p in ps}
                    ex=expected_t(L,q,side,kind)
                    check(ts==ex,f'exact whole-probability-grid fibre L={L} q={q} {side} {kind}')
                    classifications+=1
                    raycount=2*L+(len(ex)-2)*L**2
                    if L==4:
                        check(phase_count_exhaustive(ps,L,side)==raycount,f'exact phase enumeration q={q} {side} {kind}')
                    if q in (F(1,2),F(3,4)):
                        rows.append({'L':L,'q':str(q),'side':side,'grid':kind,'mixing_count':len(ts),'ray_count':raycount})
    # Independent valuation sweep does not assume the final fibre classification.
    tested=0
    for r in range(2,7):
        for n in range(1,2**r,2):
            q=F(n,2**r)
            for s in range(2,7):
                for k in range(1,2**s,2):
                    t=F(k,2**s);a=q*t+(1-q)*(1-t);b=q*t/a;c=q*(1-t)/(1-a)
                    dyadic=lambda x:x.denominator&(x.denominator-1)==0
                    check(not(dyadic(b) and dyadic(c)),'exact nonhalf dyadic mixing exclusion')
                    tested+=1
    rng=random.Random(SEED); cover_rows=[]; max_res=0.; max_step_ratio=0.
    for L in (4,8,16,64,256,4096):
        for D in sorted({1,L//4,L//2,L}):
            bound=D**-.5+4*math.sin(math.pi/(2*L)); largest=0.
            targets=[I,param(0,0,.7),param(1,0,.7)]+[random_unitary(rng) for _ in range(50)]
            for u in targets:
                v,c=quantise(u,D,L)
                err=opnorm(sub(v,u));largest=max(largest,err)
                residual=frob(sub(mm(adj(v),v),I));max_res=max(max_res,residual)
                check(residual<TOL,'numerical quantised-unitary unitarity')
                check(err<=bound+TOL,'numerical unitary covering bound')
                check((c['t']*D).denominator==1,'exact quantiser mixing certificate')
            cover_rows.append({'L':L,'D':D,'largest_sampled_operator_error':largest,'proved_upper_bound':bound})
    theta=math.pi/8;u=[[math.cos(theta)+0j,-math.sin(theta)+0j],[math.sin(theta)+0j,math.cos(theta)+0j]]
    delta=2*math.sin(math.pi/16)
    check(abs(opnorm(sub(I,u))-delta)<TOL,'numerical sharp gap attained by identity')
    for q in (F(1,8),F(1,4),F(3,4),F(7,8)):
        sq=seed_matrix(q);target=mm(u,sq)
        check(abs(chord(sq,target)-delta)<TOL,'numerical sharp state gap attained')
        for t in (0,.5,1):
            for b,g in product(range(8),repeat=2):
                v=param(t,2*math.pi*b/8,2*math.pi*g/8)
                check(chord(target,mm(v,sq))>=delta-TOL,'numerical phase sample respects sharp state gap')
    histories=[]
    for L in (16,64,256):
        for kind in ('T','K'):
            D=L if kind=='T' else L//2;bound=D**-.5+4*math.sin(math.pi/(2*L))
            v=I;actual=seed_matrix(F(1,2));ideal=actual; max_marg=0.
            for j in range(40):
                u=random_unitary(rng);side='A' if j%2==0 else 'B'
                target=mm(u,v) if side=='A' else mm(v,trans(u))
                vp,cert=quantise(target,D,L)
                local=mm(vp,adj(v)) if side=='A' else trans(mm(adj(v),vp))
                step_error=opnorm(sub(local,u));max_step_ratio=max(max_step_ratio,step_error/bound)
                check(step_error<=bound+TOL,'numerical alternating Bell local step bound')
                actual=mm(local,actual) if side=='A' else mm(actual,trans(local))
                ideal=mm(u,ideal) if side=='A' else mm(ideal,trans(u))
                v=vp
                check(frob(sub(actual,scale(v,1/math.sqrt(2))))<TOL,'numerical alternating Bell action')
                for marginal in (mm(actual,adj(actual)),mm(adj(actual),actual)):
                    res=frob(sub(marginal,scale(I,.5)));max_marg=max(max_marg,res)
                    check(res<TOL,'numerical Bell marginal preservation')
                check(cert['t'] in expected_t(L,F(1,2),'A',kind),'exact Bell endpoint probability certificate')
                check(chord(actual,ideal)<=min(math.sqrt(2),(j+1)*bound)+TOL,'numerical Bell word bound')
            histories.append({'L':L,'grid':kind,'sector':'Bell alternating A/B','steps':40,'max_marginal_residual':max_marg,'final_state_error':chord(actual,ideal)})
            q=F(3,4);D=L if kind=='T' else L//q.denominator;bound=D**-.5+4*math.sin(math.pi/(2*L))
            v=I;sq=seed_matrix(q);actual=sq;max_marg=0.
            for j in range(30):
                u=random_unitary(rng);vp,cert=quantise(mm(u,v),D,L);local=mm(vp,adj(v))
                actual=mm(actual,trans(local));v=vp
                check(opnorm(sub(local,u))<=bound+TOL,'numerical nonmaximal B-only step bound')
                check(frob(sub(actual,mm(sq,trans(v))))<TOL,'numerical nonmaximal B-only action')
                res=frob(sub(mm(actual,adj(actual)),[[.75+0j,0j],[0j,.25+0j]]));max_marg=max(max_marg,res)
                check(res<TOL,'numerical nonmaximal unchanged A marginal')
                check(cert['t'] in expected_t(L,q,'B',kind),'exact nonmaximal B-only endpoint certificate')
            histories.append({'L':L,'grid':kind,'sector':'q=3/4 B-only','steps':30,'max_marginal_residual':max_marg})
    return {'exact_fibre_classifications':classifications,'valuation_cases':tested,'fibre_rows':rows,'unitary_cover_samples':cover_rows,'sharp_gap_numeric':delta,'maximum_sampled_unitarity_residual':max_res,'maximum_history_step_to_bound_ratio':max_step_ratio,'history_rows':histories}

if __name__=='__main__':
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    out=Path(__file__).resolve().parent/'reproducibility'/('entangled_python_'+stamp+'_'+uuid.uuid4().hex[:8])
    out.mkdir(parents=True,exist_ok=False)
    record={'started_utc':datetime.now(timezone.utc).isoformat(),'script':Path(__file__).name,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python_version':sys.version,'platform':platform.platform(),'machine':platform.machine(),'numeric_tolerance':TOL,'seed':SEED,'scope':'Exact fraction/modular checks plus numerical matrix checks; no interval certification.'}
    try:
        record['results']=run();record['status']='passed'
    except Exception:
        record['status']='failed';record['exception']=traceback.format_exc()
    record['checks']=checks;record['finished_utc']=datetime.now(timezone.utc).isoformat()
    path=out/'entangled_local_run_record.json'
    path.write_text(json.dumps(record,indent=2)+'\n')
    summary=f"Verification: {record['status']}\nAssertions: {len(checks)}; passed: {sum(c['passed'] for c in checks)}\nLocal record: {path}\n"
    if record['status']=='failed':summary+=record['exception']
    (out/'entangled_local_assertion_log.txt').write_text(summary)
    print(summary)
    raise SystemExit(0 if record['status']=='passed' else 1)
