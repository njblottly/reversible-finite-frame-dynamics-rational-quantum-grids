#!/usr/bin/env python3
"""Off-diagonal marginal and finite-frame verification. Standard library only.
Run directly in PyCharm or: python3 verify_offdiagonal_frames.py
Creates fresh local JSON and text records. Exact arithmetic checks are kept
separate from numerical sampling; no interval certification is claimed.
"""
from fractions import Fraction as F
from pathlib import Path
from datetime import datetime, timezone
from itertools import product
from functools import lru_cache
import cmath
import hashlib
import json
import math
import platform
import random
import sys
import traceback
import uuid

TOL=1e-10
SEED=20260909
checks=[]
def check(test,label):
    checks.append({'label':label,'passed':bool(test)})
    if not test:raise AssertionError(label)

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


@lru_cache(None)
def sqrt_integer(n):
    """Return integer f and squarefree d with sqrt(n)=f sqrt(d)."""
    if n==0:return 0,1
    f=d=1;p=2
    while p*p<=n:
        power=0
        while n%p==0:n//=p;power+=1
        f*=p**(power//2)
        if power%2:d*=p
        p+=1
    return f,d*n

@lru_cache(None)
def row_radical(b):
    """sqrt(b(1-b)) = rational coefficient times sqrt(odd squarefree d)."""
    n,D=b.numerator,b.denominator
    f,d=sqrt_integer(n*(D-n))
    return F(f,D),d

def root_term(L,exponent,coefficient,radicand=1):
    e=exponent%L
    if e>=L//2:e-=L//2;coefficient=-coefficient
    return {} if coefficient==0 else {(radicand,e):coefficient}

def add_terms(a,b):
    c=dict(a)
    for k,v in b.items():
        c[k]=c.get(k,F(0))+v
        if c[k]==0:del c[k]
    return tuple(sorted(c.items()))

def coherence(a,b,c,u,v,L):
    x,d=row_radical(b);y,e=row_radical(c)
    return add_terms(root_term(L,-u,a*x,d),root_term(L,-v,(1-a)*y,e))

def exact_fibres(L):
    targets={}
    for r in (F(1,2),F(1,4)):
        for k in (0,1):
            targets[f'r={r},k={k}']=tuple(sorted(root_term(L,k,r/2).items()))
    targets['irrational_Schmidt_seed']=add_terms(root_term(L,0,F(1,4)),root_term(L,L//4,F(-1,4)))
    hits={name:{'T':0,'K':0} for name in targets}
    lookup={value:name for name,value in targets.items()}
    # All possible root probabilities and first conditional probabilities.
    # The second conditional probability is forced by rho_00=1/2.
    candidate_rows=0
    for m in range(1,L):
        a=F(m,L)
        for n in range(L+1):
            b=F(n,L);c=(F(1,2)-a*b)/(1-a)
            if not 0<=c<=1 or (L*c).denominator!=1:continue
            candidate_rows+=1
            p=(a*b,a*(1-b),(1-a)*c,(1-a)*(1-c))
            occupation=all((L*x).denominator==1 for x in p)
            for u,v in product(range(L),repeat=2):
                key=coherence(a,b,c,u,v,L)
                if key in lookup:
                    name=lookup[key];hits[name]['T']+=L
                    if occupation:hits[name]['K']+=L
                    check(b in (0,F(1,2),1) and c in (0,F(1,2),1),
                          'exact nonzero cyclotomic coherence row restriction')
    rows=[]
    for r in (F(1,2),F(1,4)):
        for k in (0,1):
            name=f'r={r},k={k}'
            for kind,mult in [('T',2),('K',4)]:
                expected=2*L if L%(mult*r.denominator)==0 else 0
                check(hits[name][kind]==expected,f'exact full row/phase enumeration {name} L={L} {kind}')
                rows.append({'L':L,'family':name,'grid':kind,'rays':hits[name][kind]})
    for kind in ('T','K'):
        check(hits['irrational_Schmidt_seed'][kind]==2*L,'exact irrational-Schmidt seed fibre count')
        rows.append({'L':L,'family':'irrational_Schmidt_seed','grid':kind,'rays':hits['irrational_Schmidt_seed'][kind]})
    return rows,candidate_rows

def prepare(phi,L,kind):
    probs=[abs(x)**2 for row in phi for x in row]
    a=probs[0]+probs[1];b=probs[0]/a if a else 0.;c=probs[2]/(1-a) if a!=1 else 0.
    if kind=='T':A=B=L
    else:A=2**((L.bit_length()-1)//2);B=L//A
    aa,bb,cc=F(round(A*a),A),F(round(B*b),B),F(round(B*c),B)
    pp=[aa*bb,aa*(1-bb),(1-aa)*cc,(1-aa)*(1-cc)]
    values=[x for row in phi for x in row];ref=next(x for x in values if abs(x)>0)
    phases=[round(cmath.phase(x/ref)*L/(2*math.pi))%L if abs(x)>0 else 0 for x in values]
    vals=[math.sqrt(float(p))*cmath.exp(2j*math.pi*h/L) for p,h in zip(pp,phases)]
    seed=[vals[:2],vals[2:]]
    check(all((L*x).denominator==1 for x in (aa,bb,cc)),'exact preparation conditional certificate')
    if kind=='K':check(all((L*x).denominator==1 for x in pp),'exact preparation occupation certificate')
    eta=math.sqrt(1/A+1/B)+2*math.sin(math.pi/(2*L))
    check(chord(phi,seed)<=eta+TOL,'numerical preparation covering bound')
    return seed,eta

def run():
    rows=[];candidate_rows=0
    for L in (4,8,16):
        rr,n=exact_fibres(L);rows+=rr;candidate_rows+=n
    for M in range(3,13):
        L=2**M
        powers={pow(5,j,L) for j in range(2**(M-2))}
        check(len(powers)==2**(M-2),'exact order of 5 modulo dyadic L')
        check(powers|{(-x)%L for x in powers}==set(range(1,L,2)),
              'exact generators of dyadic cyclotomic Galois group')
    for s in range(2,11):
        for n in range(1,2**s,2):
            N=n*(2**s-n)
            check(N%4==3 and math.isqrt(N)**2!=N,'exact nonhalf conditional radicand obstruction')
    h=scale([[1+0j,1+0j],[1+0j,-1+0j]],1/math.sqrt(2));z=[[1+0j,0j],[0j,-1+0j]]
    sq=seed_matrix(F(3,4));omega=mm(sq,trans(h));target=mm(h,omega)
    gap=math.sqrt(2-math.sqrt(2))
    check(abs(opnorm(sub(z,h))-gap)<TOL,'numerical offdiagonal sharp operator gap attained')
    check(abs(chord(mm(z,omega),target)-gap)<TOL,'numerical offdiagonal sharp state gap attained')
    for L in (4,8,16):
        for phase in range(L):
            diag=[[1+0j,0j],[0j,cmath.exp(2j*math.pi*phase/L)]]
            for v in (diag,mm(diag,[[0j,1+0j],[1+0j,0j]])):
                check(chord(mm(v,omega),target)>=gap-TOL,'numerical monomial phases respect offdiagonal state gap')
    nu=scale([[1+0j,1+0j],[1+0j,1j]],.5);rho=mm(nu,adj(nu))
    check(abs(rho[0][1]-(1-1j)/4)<TOL,'numerical irrational-Schmidt marginal coherence')
    check(abs((rho[0][0]*rho[1][1]-rho[0][1]*rho[1][0])-F(1,8))<TOL,
          'numerical irrational-Schmidt determinant')
    idealp=[(2+math.sqrt(3))/8,(2-math.sqrt(3))/8,(2-math.sqrt(3))/8,(2+math.sqrt(3))/8]
    check(max(abs(abs(x)**2-p) for x,p in zip([x for row in target for x in row],idealp))<TOL,
          'numerical atlas output outside the fixed grid witness')
    # Tradeoff: choose changed preparations in the same exact marginal fibre;
    # selected output is C times a fixed-grid fibre state. Any chart C is allowed here.
    rng=random.Random(SEED)
    for _ in range(30):
        W=random_unitary(rng);C=random_unitary(rng);T=z
        V=mm(mm(C,T),adj(W));psip=mm(W,omega)
        total=chord(omega,psip)+opnorm(sub(V,h))+opnorm(sub(C,I))
        check(total>=gap-TOL,'numerical fixed-marginal tolerance tradeoff (unminimised operator distances)')
    history_rows=[];max_marg=0.;max_unitarity=0.
    for L in (16,64,256):
        D=L;eps=D**-.5+4*math.sin(math.pi/(2*L))
        for kind in ('T','K'):
            vals=[complex(rng.uniform(-1,1),rng.uniform(-1,1)) for _ in range(4)]
            norm=math.sqrt(sum(abs(x)**2 for x in vals));vals=[x/norm for x in vals]
            nominal=[vals[:2],vals[2:]];rounded,eta=prepare(nominal,L,kind)
            for label,seed,ideal,initial_error in [('irrational Schmidt',nu,nu,0.),('selected preparation',rounded,nominal,eta)]:
                A=I;B=I;actual=seed;worst=0.;worststep=0.
                for j in range(30):
                    U=random_unitary(rng);side='A' if j%2==0 else 'B'
                    old_remote=mm(adj(actual),actual) if side=='A' else mm(actual,adj(actual))
                    if side=='A':
                        Ap,cert=quantise(mm(U,A),D,L);V=mm(Ap,adj(A));A=Ap
                        actual=mm(V,actual);ideal=mm(U,ideal)
                    else:
                        Bp,cert=quantise(mm(U,B),D,L);V=mm(Bp,adj(B));B=Bp
                        actual=mm(actual,trans(V));ideal=mm(ideal,trans(U))
                    err=opnorm(sub(V,U));worststep=max(worststep,err)
                    check(err<=eps+TOL,'numerical finite-frame step bound')
                    check((D*cert['t']).denominator==1,'exact finite-frame probability label')
                    for Fm in (A,B,V):
                        res=frob(sub(mm(adj(Fm),Fm),I));max_unitarity=max(max_unitarity,res)
                        check(res<TOL,'numerical frame and selected-gate unitarity')
                    check(frob(sub(actual,mm(mm(A,seed),trans(B))))<TOL,'numerical atlas endpoint label reconstruction')
                    remote=mm(adj(actual),actual) if side=='A' else mm(actual,adj(actual))
                    res=frob(sub(remote,old_remote));max_marg=max(max_marg,res);worst=max(worst,res)
                    check(res<TOL,'numerical remote marginal unchanged at local step')
                    check(chord(actual,ideal)<=min(math.sqrt(2),initial_error+(j+1)*eps)+TOL,'numerical atlas history/preparation bound')
                history_rows.append({'L':L,'grid_seed':kind,'seed':label,'steps':30,'maximum_step_error':worststep,'step_upper_bound':eps,'maximum_remote_marginal_residual':worst,'final_state_error':chord(actual,ideal)})
            # Check reordered independent local controls produce the same frame-labelled endpoint.
            UA=random_unitary(rng);UB=random_unitary(rng)
            A1,_=quantise(mm(UA,A),D,L);B1,_=quantise(mm(UB,B),D,L)
            VA=mm(A1,adj(A));VB=mm(B1,adj(B));state=mm(mm(A,nu),trans(B))
            ab=mm(mm(VA,state),trans(VB));ba=mm(VA,mm(state,trans(VB)))
            check(frob(sub(ab,ba))<TOL,'numerical commuting distinct-subsystem frame updates')
    return {'fibre_rows':rows,'candidate_probability_rows':candidate_rows,
      'sharp_gap_numeric':gap,'maximum_unitarity_residual':max_unitarity,
      'maximum_remote_marginal_residual':max_marg,'history_rows':history_rows}
if __name__=='__main__':
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    out=Path(__file__).resolve().parent/'reproducibility'/('offdiagonal_python_'+stamp+'_'+uuid.uuid4().hex[:8])
    out.mkdir(parents=True,exist_ok=False)
    record={'started_utc':datetime.now(timezone.utc).isoformat(),'script':Path(__file__).name,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python_version':sys.version,'platform':platform.platform(),'machine':platform.machine(),'numeric_tolerance':TOL,'seed':SEED,'scope':'Exact fraction/modular checks plus numerical matrix checks; no interval certification.'}
    try:
        record['results']=run();record['status']='passed'
    except Exception:
        record['status']='failed';record['exception']=traceback.format_exc()
    record['checks']=checks;record['finished_utc']=datetime.now(timezone.utc).isoformat()
    path=out/'offdiagonal_frames_run_record.json'
    path.write_text(json.dumps(record,indent=2)+'\n')
    summary=f"Verification: {record['status']}\nAssertions: {len(checks)}; passed: {sum(c['passed'] for c in checks)}\nLocal record: {path}\n"
    if record['status']=='failed':summary+=record['exception']
    (out/'offdiagonal_frames_assertion_log.txt').write_text(summary)
    print(summary)
    raise SystemExit(0 if record['status']=='passed' else 1)
