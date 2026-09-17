#!/usr/bin/env python3
"""Reversible finite-frame checks, version 0.1. Python 3.10+, standard library.
Exact finite-map/cell-volume checks and separate numerical matching checks.
Each run saves new JSON and text records; no interval certification claimed.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from datetime import datetime, timezone
import cmath, math, random, json, hashlib, platform, sys, traceback, uuid
TOL=1e-9
SEED=20260911
checks=[]
def check(test,label):
    checks.append({'label':label,'passed':bool(test)})
    if not test: raise AssertionError(label)
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



def catalogue(D,L):
    # Each ray has exactly one label, including collapsed endpoint phases.
    out=[param(0,0,2*math.pi*k/L) for k in range(L)]
    out += [param(1,0,2*math.pi*k/L) for k in range(L)]
    out += [param(n/D,2*math.pi*b/L,2*math.pi*g/L)
            for n in range(1,D) for b in range(L) for g in range(L)]
    return out

def pdist(a,b):
    trace=sum(x.conjugate()*y for ar,br in zip(a,b) for x,y in zip(ar,br))
    if abs(trace)<1e-12: return math.sqrt(2)
    return opnorm(sub(a,scale(b,trace.conjugate()/abs(trace))))

def matching(cost,threshold):
    n=len(cost); right=[-1]*n
    edges=[[j for j in range(n) if cost[i][j]<=threshold] for i in range(n)]
    def augment(i,seen):
        for j in edges[i]:
            if seen[j]:continue
            seen[j]=True
            if right[j]<0 or augment(right[j],seen):
                right[j]=i;return True
        return False
    for i in range(n):
        if not augment(i,[False]*n):return None
    p=[-1]*n
    for j,i in enumerate(right):p[i]=j
    return p

def bottleneck(cost):
    vals=sorted(set(x for row in cost for x in row));lo=-1;hi=len(vals)-1
    while hi-lo>1:
        mid=(hi+lo)//2
        if matching(cost,vals[mid]) is None:lo=mid
        else:hi=mid
    return matching(cost,vals[hi]),vals[hi],(vals[lo] if lo>=0 else None)

def inverse(p):
    q=[0]*len(p)
    for i,j in enumerate(p):q[j]=i
    return q

def clean_extension(q):
    n=len(q);fibres=[[x for x in range(n) if q[x]==y] for y in range(n)]
    m=max(map(len,fibres));p={}
    for group in fibres:
        for rank,x in enumerate(group):p[x*m]=q[x]*m+rank
    unused_in=[z for z in range(n*m) if z not in p]
    used=set(p.values());unused_out=[z for z in range(n*m) if z not in used]
    p.update(zip(unused_in,unused_out))
    return m,[p[z] for z in range(n*m)]

def run():
    rng=random.Random(SEED);rows=[];history=[];maxres=0.
    map_count=0
    for n in range(1,5):
        for q in product(range(n),repeat=n):
            m,p=clean_extension(q);map_count+=1
            check(sorted(p)==list(range(n*m)),'exact clean-ancilla extension is a permutation')
            check(all(p[x*m]//m==q[x] for x in range(n)),'exact clean-input rounding reproduced')
            check(m==max(q.count(y) for y in range(n)),'exact largest-fibre auxiliary alphabet')
    for L in (4,8,16,32,64,256):
        for D in (2,L):
            N=2*L+(D-1)*L*L;h=F(D-1)+F(2,L);c=1/h;a=c/L
            check(2*a+(D-1)*c==1,'exact Haar slab volumes sum to one')
            check(c/(L*L)==F(1,N) and a/L==F(1,N),'exact all frame cells have equal Haar volume')
            check(all(a+(n-1)*c<=F(n,D)<=a+n*c for n in range(1,D)),
                  'exact interior frame lies in its assigned slab')
    gates=[('Hadamard',param(.5,0,0)),('real rotation',
        [[math.cos(math.pi/7),-math.sin(math.pi/7)],[math.sin(math.pi/7),math.cos(math.pi/7)]]),
        ('complex gate',param(1/3,math.pi/5,math.pi/7))]
    nu=[[.5+0j,.5+0j],[.5+0j,.5j]]
    for D,L in ((2,4),(4,4),(4,8),(8,8)):
        cat=catalogue(D,L);N=len(cat);B=2*(math.sqrt(2/(D-1+2/L))+4*math.sin(math.pi/(2*L)))
        check(N==2*L+(D-1)*L*L,'exact catalogue size')
        # Numerical samples of the analytic cell-radius estimate, including poles.
        c=1/(D-1+2/L);a=c/L;r=B/2
        for _ in range(60):
            n=rng.randrange(D+1);b=rng.randrange(L);g=rng.randrange(L)
            if n in (0,D):
                t=rng.random()*a if n==0 else 1-rng.random()*a
                beta=rng.uniform(-math.pi,math.pi);theta=2*math.pi*g/L+rng.uniform(-math.pi/L,math.pi/L)
                gamma=theta+beta if n==0 else theta-beta
                centre=param(n/D,0,2*math.pi*g/L)
            else:
                t=a+(n-1+rng.random())*c
                beta=2*math.pi*b/L+rng.uniform(-math.pi/L,math.pi/L)
                gamma=2*math.pi*g/L+rng.uniform(-math.pi/L,math.pi/L)
                centre=param(n/D,2*math.pi*b/L,2*math.pi*g/L)
            check(pdist(centre,param(t,beta,gamma))<=r+TOL,'numerical cell-radius sample')
        tables=[]
        for name,U in gates:
            targets=[mm(U,A) for A in cat]
            cost=[[pdist(A,T) for A in cat] for T in targets]
            p,b,prev=bottleneck(cost);pinv=inverse(p)
            check(sorted(p)==list(range(N)),'exact selected label map is a permutation')
            check(all(pinv[p[i]]==i for i in range(N)),'exact label inverse')
            check(max(cost[i][p[i]] for i in range(N))<=b+TOL,'numerical bottleneck edge bound')
            check(prev is None or matching(cost,prev) is None,'numerical cost graph has no matching at preceding threshold')
            check(b<=B+TOL,'numerical matching obeys analytic uniform bound')
            worstinverse=max(pdist(cat[pinv[j]],mm(adj(U),cat[j])) for j in range(N))
            check(worstinverse<=b+TOL,'numerical inverse command has same error bound')
            rows.append({'D':D,'L':L,'N':N,'gate':name,'bottleneck_numeric':b,'proved_upper_bound':B})
            tables.append((U,p,pinv,b))
        ai=L+L//2;bi=ai # t=1, gamma=pi represents the identity.
        start=(ai,bi);actual=nu;ideal=nu;budget=0;word=[]
        for step in range(24):
            gi=step%3;side=step%2;U,p,pinv,b=tables[gi]
            previous=mm(adj(actual),actual) if side==0 else mm(actual,adj(actual))
            if side==0:
                nxt=p[ai];V=mm(cat[nxt],adj(cat[ai]));ai=nxt
                actual=mm(V,actual);ideal=mm(U,ideal)
            else:
                nxt=p[bi];V=mm(cat[nxt],adj(cat[bi]));bi=nxt
                actual=mm(actual,trans(V));ideal=mm(ideal,trans(U))
            word.append((side,gi));budget+=b
            res=frob(sub(mm(adj(V),V),I));maxres=max(maxres,res)
            check(res<TOL,'numerical selected local matrix unitary')
            remote=mm(adj(actual),actual) if side==0 else mm(actual,adj(actual))
            check(frob(sub(remote,previous))<TOL,'numerical remote marginal preserved')
            check(frob(sub(actual,mm(mm(cat[ai],nu),trans(cat[bi]))))<TOL,'numerical labelled endpoint')
            check(chord(actual,ideal)<=min(math.sqrt(2),budget)+TOL,'numerical reversible word error bound')
        final_error=chord(actual,ideal)
        for side,gi in reversed(word):
            _,_,pinv,_=tables[gi]
            if side==0:
                nxt=pinv[ai];V=mm(cat[nxt],adj(cat[ai]));ai=nxt;actual=mm(V,actual)
            else:
                nxt=pinv[bi];V=mm(cat[nxt],adj(cat[bi]));bi=nxt;actual=mm(actual,trans(V))
        check((ai,bi)==start,'exact full word reversal restores frame labels')
        check(chord(actual,nu)<TOL,'numerical full word reversal restores physical state')
        p=tables[0][1];q=tables[1][1]
        fa=lambda pair:(p[pair[0]],pair[1])
        fb=lambda pair:(pair[0],q[pair[1]])
        check(all(fa(fb((i,j)))==fb(fa((i,j))) for i in range(N) for j in range(N)),
              'exact independent coordinate updates commute on every frame pair')
        history.append({'D':D,'L':L,'forward_steps':24,'final_forward_error':final_error,'return_residual':chord(actual,nu)})
    return {'finite_maps_exhausted':map_count,'matching_rows':rows,'histories':history,'maximum_unitarity_residual':maxres}

if __name__=='__main__':
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    out=Path(__file__).resolve().parent/'reproducibility'/('reversible_python_'+stamp+'_'+uuid.uuid4().hex[:8]);out.mkdir(parents=True)
    r={'started_utc':datetime.now(timezone.utc).isoformat(),'script':Path(__file__).name,'script_version':'0.1',
       'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python_version':sys.version,
       'platform':platform.platform(),'seed':SEED,'numeric_tolerance':TOL,
       'scope':'Exact combinatorial checks and separately numerical matrix/matching checks; not interval certification.'}
    try:r['results']=run();r['status']='passed'
    except Exception:r['status']='failed';r['exception']=traceback.format_exc()
    r['checks']=checks;r['finished_utc']=datetime.now(timezone.utc).isoformat()
    path=out/'reversible_frames_run_record.json';path.write_text(json.dumps(r,indent=2)+'\n')
    text=f"Verification: {r['status']}\nAssertions: {len(checks)}; passed: {sum(c['passed'] for c in checks)}\nLocal record: {path}\n"
    if r['status']=='failed':text+=r['exception']
    (out/'reversible_frames_assertion_log.txt').write_text(text);print(text)
    raise SystemExit(0 if r['status']=='passed' else 1)
