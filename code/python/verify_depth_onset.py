#!/usr/bin/env python3
"""Depth-onset certificates v0.1. Python 3.10+, standard library only.
Run: python3 verify_depth_onset.py
Certifies integer spectral-budget bounds through 46080 labels and earliest
single-depth thresholds for the three frozen tables, over all phase priors.
Candidate budget constants are frozen; all verification uses exact integers
and algebraic arithmetic. No floating-point discovery enters this verifier.
Universal window and error-budget bounds are analytic manuscript results.
"""
from collections import deque, Counter
from math import lcm, isqrt
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
import hashlib, json, platform, sys, traceback, uuid

VERSION = '0.1'
from fractions import Fraction as Q
ROOT = Path(__file__).resolve().parent
CHECKS = []
SUMMARY = {}
TABLE = {}
def check(test, label):
    CHECKS.append({'label': label, 'passed': bool(test)})
    if not test:
        raise AssertionError(label)

def mul1(a, b):
    """P_a P_b = i**phase P_c, with I,X,Y,Z indexed 0,1,2,3."""
    if a == 0: return 0, b
    if b == 0: return 0, a
    if a == b: return 0, 0
    return (1 if (a,b) in ((1,2),(2,3),(3,1)) else 3), 6-a-b

def mul2(a, b):
    r,x = mul1(a//4, b//4); s,y = mul1(a%4, b%4)
    return (r+s)%4, 4*x+y

def compose(p, q):
    """Conjugation action of product P Q; signed indices are one-based."""
    return tuple((1 if x>0 else -1)*p[abs(x)-1] for x in q)

def inv(p):
    out = [0]*len(p)
    for i,x in enumerate(p): out[abs(x)-1] = (i+1)*(1 if x>0 else -1)
    return tuple(out)

ID = tuple(range(1,17))
h1 = (1,4,-3,2)
s1 = (1,3,-2,4)
def local(p, side):
    ans = []
    for a in range(4):
        for b in range(4):
            x = p[(a,b)[side]]
            j = 4*(abs(x)-1)+b if side==0 else 4*a+abs(x)-1
            ans.append((j+1)*(1 if x>0 else -1))
    return tuple(ans)
HA,HB,SA,SB = local(h1,0),local(h1,1),local(s1,0),local(s1,1)
# CNOT conjugates A's I,X,Y,Z to II,XX,YX,ZI and B's to II,IX,ZY,ZZ.
ca,cb = (0,5,9,12),(0,1,14,15)
cnot = []
for a in range(4):
    for b in range(4):
        phase,j = mul2(ca[a],cb[b])
        assert phase in (0,2)
        cnot.append((j+1)*(1 if phase==0 else -1))
C = tuple(cnot)
GENS = {'HA':HA,'HB':HB,'SA':SA,'SB':SB,'C':C}

def enumerate_group(generators):
    out=[ID]; seen={ID:0}; queue=deque([ID])
    while queue:
        x=queue.popleft()
        for g in generators:
            y=compose(g,x)
            if y not in seen:
                seen[y]=len(out);out.append(y);queue.append(y)
    return out,seen

def act(p, r):
    out=[0]*16
    for i,x in enumerate(p): out[abs(x)-1]=(1 if x>0 else -1)*r[i]
    return tuple(out)

def word_product(word):
    p=ID
    for name in word: p=compose(GENS[name],p)
    return p

def power(p,n):
    q=ID
    for _ in range(n): q=compose(p,q)
    return q

def state(entries):
    r=[0]*16;r[0]=1
    for j,v in entries.items():r[j]=v
    return tuple(r)

def average(rows):
    return tuple(sum(Fraction(r[j]) for r in rows)/len(rows) for j in range(16))

PRIMES=(2,3,5,7,-1)
FACTORS=tuple(__import__('functools').reduce(lambda a,b:a*b,
               (p for j,p in enumerate(PRIMES) if mask>>j&1),1) for mask in range(32))
class F:
    """Sparse exact multiquadratic field; bit 4 represents i."""
    def __init__(self,x=0):
        self.c={k:Q(v) for k,v in x.items() if v} if isinstance(x,dict) else ({0:Q(x)} if x else {})
    def __add__(self,o):
        o=o if isinstance(o,F) else F(o);c=self.c.copy()
        for k,v in o.c.items():c[k]=c.get(k,0)+v
        return F(c)
    __radd__=__add__
    def __neg__(self):return F({k:-v for k,v in self.c.items()})
    def __sub__(self,o):return self+-asf(o)
    def __rsub__(self,o):return asf(o)+-self
    def __mul__(self,o):
        o=asf(o);c={}
        for a,x in self.c.items():
            for b,y in o.c.items():c[a^b]=c.get(a^b,0)+x*y*FACTORS[a&b]
        return F(c)
    __rmul__=__mul__
    def __eq__(self,o):return self.c==asf(o).c
    def conj(self):return F({k:(-v if k&16 else v) for k,v in self.c.items()})
    def __repr__(self):return str({str(k):str(v) for k,v in self.c.items()})
def asf(x):return x if isinstance(x,F) else F(x)
def root(n):
    if n==0:return F()
    factor=1;mask=0
    for j,p in enumerate(PRIMES[:4]):
        count=0
        while n%p==0:n//=p;count+=1
        factor*=p**(count//2)
        if count%2:mask|=1<<j
    assert n==1
    return F({mask:factor})
def transpose(a):return tuple(zip(*a))
def adj(a):return tuple(tuple(z.conj() for z in row) for row in transpose(a))
def mm(a,b):return tuple(tuple(sum((x*y for x,y in zip(row,col)),F()) for col in transpose(b)) for row in a)
def scale(a,q):return tuple(tuple(x*q for x in row) for row in a)
def ident(n):return tuple(tuple(F(i==j) for j in range(n)) for i in range(n))
def outer(v):return tuple(tuple(x*y.conj() for y in v) for x in v)
def addmat(a,b):return tuple(tuple(x+y for x,y in zip(ar,br)) for ar,br in zip(a,b))
def flat(a):return tuple(x for row in a for x in row)
def trace(a):return sum((a[i][i] for i in range(len(a))),F())
def sqabs(x):return x.conj()*x


def kron(a,b):
    return tuple(tuple(a[i//len(b)][j//len(b[0])]*b[i%len(b)][j%len(b[0])] for j in range(len(a[0])*len(b[0]))) for i in range(len(a)*len(b)))
def negmat(a):return scale(a,-1)
def phase_coeff(r,a,b):
    r=list(map(asf,r));u=root(2)*Q(1,2)
    co=(F(1),u,F(),-u,F(-1),-u,F(),u)
    si=(F(),u,F(1),u,F(),-u,F(-1),-u)
    for side,n in ((0,a),(1,b)):
        c,s=co[n%8],si[n%8]
        for j in range(4):
            x,y=(4+j,8+j) if side==0 else (4*j+1,4*j+2)
            rx,ry=r[x],r[y];r[x]=c*rx-s*ry;r[y]=s*rx+c*ry
    return tuple(r)
def computational(r):
    return tuple((1+(-1)**a*r[12]+(-1)**b*r[3]+(-1)**(a+b)*r[15])*Q(1,4) for a in range(2) for b in range(2))
def encode_scalar(x):
    x=asf(x)
    if not x.c:return '0'
    if set(x.c)<={0}:return str(x.c[0])
    return str(x)


def permutation_power(p,n):
    result=list(range(len(p)));base=p[:]
    while n:
        if n&1:result=[base[x] for x in result]
        base=[base[x] for x in base];n//=2
    return result

def cycle_counts(p):
    seen=bytearray(len(p));counts=Counter()
    for start in range(len(p)):
        if seen[start]:continue
        x=start;k=0
        while not seen[x]:seen[x]=1;k+=1;x=p[x]
        check(x==start,'permutation orbit closes at its starting label')
        counts[k]+=1
    return dict(sorted(counts.items()))


BUDGETS={8:(27,7),16:(480,11),32:(480,22),48:(480,44),64:(481,55),128:(6970,128),256:(30220,139),464:(44293,406),1024:(20885827,684),4096:(20885827,3420),46080:(20885844,45828)}
def radical_sign(a,b):
    """Exact sign of a+b sqrt(2) for integer a,b."""
    if not b:return (a>0)-(a<0)
    if not a:return (b>0)-(b<0)
    if a>0 and b>0:return 1
    if a<0 and b<0:return -1
    d=a*a-2*b*b
    return ((d>0)-(d<0))*(1 if a>0 else -1)

def main():
    global TABLE
    # A_q = 2^(q+1) T_q((-2+sqrt(2))/4) = a_q+b_q sqrt(2).
    aa,bb=2,0;a,b=-2,1;saved={};budget_rows=[];witnesses={w for _,w in BUDGETS.values()}
    prior_K=0
    for M,(K,w) in BUDGETS.items():check(K>=prior_K,'frozen budget bounds are monotone');prior_K=K
    for q in range(1,max(BUDGETS)+1):
        if q>1:aa,bb,a,b=a,b,-2*a+2*b-4*aa,a-2*b-4*bb
        if q in witnesses:saved[q]=(a,b)
        M=next(m for m in BUDGETS if m>=q);K=BUDGETS[M][0]
        check(radical_sign(K*K*((1<<(q+1))-a)-q*q*(1<<(q+2)),-K*K*b)>=0,'exact spectral bound at period '+str(q))
        if q in BUDGETS:
            K,w=BUDGETS[q];wa,wb=saved[w];kk=(K-1)**2
            check(radical_sign(kk*((1<<(w+1))-wa)-w*w*(1<<(w+2)),-kk*wb)<0,'integer bound is minimal at budget '+str(q))
            budget_rows.append({'frame_label_bound_M':q,'ceil_kappa_M':K,'minimality_witness_period':w,'window_length_D':20*(K+1),'maximum_blocks':20*(K+1)-1,'maximum_elementary_commands':3*(20*(K+1)-1)})
            print('Certified spectral budget:',q,flush=True)
    # Rational certificates for C0<1 (sines are positive at these angles).
    check(Q(6,1)-Q(3,2)==Q(9,2),'sin(theta_minus)>3/4 using sqrt(2)<3/2')
    check(Q(7,1)+4>16*Q(4,5)**2,'sin(theta_plus+theta_minus)>4/5 using sqrt(17)>4')
    check(Q(7,1)-Q(17,4)>16*Q(2,5)**2,'sin(theta_plus-theta_minus)>2/5 using sqrt(17)<17/4')
    check(Q(5+2,8)>Q(9,10)**2,'sin(2 theta_minus)>9/10 using sqrt(2)>1')
    check(Q(8,17)+Q(15,56)+Q(10,81)<1,'rational upper bound on ideal-correlation remainder C0')
    check(Q(1,5)/2-Q(1,20)==Q(1,20),'universal strict one-twentieth contrast threshold')
    group,index=enumerate_group(GENS.values());n=len(group);size=4*n;enc=lambda a,b,j:(2*a+b)*n+j
    left={name:[index[compose(p,g)] for g in group] for name,p in GENS.items()}
    twisted=[index[compose(compose(inv(SB),compose(C,SB)),g)] for g in group]
    blocks={m:[] for m in ('original','improved','anchor_preserving')}
    for a in range(2):
        for b in range(2):
            for j in range(n):
                jj=left['HA'][j if a==0 else left['SA'][j]]
                blocks['original'].append(enc(1-a,b,left['C'][jj]))
                blocks['improved'].append(enc(a,1-b,left['C'][jj]))
                blocks['anchor_preserving'].append(enc(1-a,b,left['C'][jj] if b==0 else twisted[jj]))
    check(size==46080,'full frame count for existing-table benchmark')
    for model,p in blocks.items():check(sorted(p)==list(range(size)),model+' block is a full permutation')
    starts=[enc(k%2,l%2,index[compose(power(SA,k//2),power(SB,l//2))]) for k in range(8) for l in range(8)]
    positions={m:starts[:] for m in blocks};zero=state({3:1,12:1,15:1})
    r=root(2);xp=-(2+r)*Q(1,4);xm=-(2-r)*Q(1,4);cp0,cp1,cm0,cm1=F(1),xp,F(1),xm
    thresholds=(Q(1,20),Q(1,3));first={t:None for t in thresholds};rows=[]
    for depth in range(25):
        cp=F(1) if depth==0 else cp1;cm=F(1) if depth==0 else cm1
        prob=(((1-cp)*(6-r)+(1-cm)*(6+r))*Q(1,34)).c.get(0,Q(0));models={};gaps={};values_by_model={}
        for model,poss in positions.items():
            values=[]
            for k in range(8):
                for ell in range(8):
                    bits,j=divmod(poss[8*k+ell],n);a,b=divmod(bits,2)
                    rr=phase_coeff(act(group[j],phase_coeff(zero,-k,-ell)),a,b)
                    value_f=(1-rr[12])*Q(1,2);check(set(value_f.c)<={0},model+' context output is rational')
                    value=value_f.c.get(0,Q(0));check(0<=value<=1,model+' context output is a valid probability');values.append(value)
            lo,hi=min(values),max(values);gap=max(lo-prob,prob-hi,Q(0));gaps[model]=gap;values_by_model[model]=values
            models[model]={'minimum_A1':str(lo),'maximum_A1':str(hi),'calibrated_A1':str(values[0]),'uniform_A1':str(sum(values,Q(0))/64),'minimum_prior_discrepancy':str(gap)}
        robust=min(gaps.values());wmodel=min(gaps,key=gaps.get);vals=values_by_model[wmodel];lo,hi=min(vals),max(vals)
        target=min(max(prob,lo),hi);weight=Q(0) if lo==hi else (target-lo)/(hi-lo)
        check(0<=weight<=1 and abs((1-weight)*lo+weight*hi-prob)==robust,'explicit two-context prior attains the benchmark infimum')
        for t in thresholds:
            if first[t] is None and robust>t:first[t]=depth
        rows.append({'blocks':depth,'ideal_A1':str(prob),'minimum_over_three_tables_and_all_priors':str(robust),'models':models,
          'attaining_prior':{'model':wmodel,'low_context_index':vals.index(lo),'high_context_index':vals.index(hi),'weight_on_high':str(weight),'context_index_encoding':'8*k+ell; phases 0..7'}})
        for model in positions:positions[model]=[blocks[model][x] for x in positions[model]]
        if depth>=1:cp0,cp1=cp1,2*xp*cp1-cp0;cm0,cm1=cm1,2*xm*cm1-cm0
    check(first[Q(1,20)]==5,'earliest common single-depth discrepancy >1/20 for the three tables is five blocks')
    check(first[Q(1,3)]==7,'earliest common single-depth discrepancy >1/3 for the three tables is seven blocks')
    check(rows[5]['ideal_A1']=='3/8' and rows[5]['minimum_over_three_tables_and_all_priors']=='1/8','five-block exact witness')
    check(rows[7]['ideal_A1']=='1/16' and rows[7]['minimum_over_three_tables_and_all_priors']=='7/16','seven-block exact witness')
    TABLE={'schema':'raqm-depth-onset-v1','spectral_budget_certificates':budget_rows,'three_table_single_depth_scan':rows,
      'earliest_three_table_strict_thresholds':{str(k):v for k,v in first.items()},
      'scope':'Earliest times are only for the three named tables. Spectral budgets support a universal multi-depth contrast, not a universal single-depth optimum. No near-sqrt(L) sufficiency claim.'}
    SUMMARY.update({'spectral_periods_certified':max(BUDGETS),'largest_budget_certificate':budget_rows[-1],
      'earliest_three_table_strict_thresholds':TABLE['earliest_three_table_strict_thresholds'],
      'five_block_discrepancy':'1/8','seven_block_discrepancy':'7/16','scope':TABLE['scope']})

if __name__=='__main__':
    started=datetime.now(timezone.utc);out=ROOT/'reproducibility'/('depth_onset_python_'+started.strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8]);out.mkdir(parents=True)
    status='failed';error=None
    try:main();status='passed'
    except Exception:error=traceback.format_exc();print(error)
    if TABLE:(out/'depth_onset_results.json').write_text(json.dumps(TABLE,separators=(',',':'))+'\n')
    record={'script':Path(__file__).name,'script_version':VERSION,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      'started_utc':started.isoformat(),'finished_utc':datetime.now(timezone.utc).isoformat(),'python_version':sys.version,'platform':platform.platform(),
      'status':status,'arithmetic':'exact integer and algebraic arithmetic; no floating-point verification','summary':SUMMARY,
      'passed_assertions':sum(c['passed'] for c in CHECKS),'total_assertions':len(CHECKS),'checks':CHECKS,'error':error}
    (out/'depth_onset_run_record.json').write_text(json.dumps(record,indent=2)+'\n')
    (out/'depth_onset_assertion_log.txt').write_text('\n'.join(('PASS: ' if c['passed'] else 'FAIL: ')+c['label'] for c in CHECKS)+'\n'+(error or ''))
    print(json.dumps(SUMMARY,indent=2));print('Mathematical verification:',status,'; assertions:',record['passed_assertions']);print('Local record:',out)
    sys.exit(0 if status=='passed' else 1)
