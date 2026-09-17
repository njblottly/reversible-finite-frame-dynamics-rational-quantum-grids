#!/usr/bin/env python3
"""Repeated mixed-control recurrence verifier v0.1. Python 3.10+, standard library.
Run: python3 verify_mixed_depth.py
Exact finite checks for three 46080-frame models, ideal entangling probabilities,
and algebraic ingredients of the universal recurrence theorem.
The arbitrary-resolution and infinite-time claims are analytic proofs in the
manuscript, not numerical extrapolations. No large refined table is built here.
Fresh JSON, assertion log and result table are saved beside this script.
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

def main():
    global TABLE
    group,index=enumerate_group(GENS.values());n=len(group);size=4*n
    check(n==11520,'Clifford base has 11520 exact signed Pauli actions')
    enc=lambda a,b,g:(2*a+b)*n+g
    def dec(x):q,g=divmod(x,n);return q//2,q%2,g
    left={name:[index[compose(p,g)] for g in group] for name,p in GENS.items()}
    twisted=[index[compose(compose(inv(SB),compose(C,SB)),g)] for g in group]
    names=('HA','HB','TA','TB','C');old={name:[] for name in names};new={name:[] for name in names};anchor_c=[]
    for a in range(2):
        for b in range(2):
            for g in range(n):
                old['HA'].append(enc(a,b,left['HA'][g]));new['HA'].append(enc(1-a,b,left['HA'][g]))
                old['HB'].append(enc(a,b,left['HB'][g]));new['HB'].append(enc(a,1-b,left['HB'][g]))
                old['C'].append(enc(a,b,left['C'][g]));new['C'].append(enc(a,1-b,left['C'][g]))
                anchor_c.append(enc(a,b,left['C'][g] if b==0 else twisted[g]))
                old['TA'].append(enc(1-a,b,g if a==0 else left['SA'][g]))
                old['TB'].append(enc(a,1-b,g if b==0 else left['SB'][g]))
    new['TA']=old['TA'];new['TB']=old['TB'];anchored={**old,'C':anchor_c}
    variants={'original':old,'improved':new,'anchor_preserving':anchored}
    blocks={};cycles={};orders={};hashes={}
    expected_cycles={'original':{16:2880},'improved':{8:2880,12:1920},'anchor_preserving':{16:2880}}
    for model,tables in variants.items():
        hashes[model]={}
        for name,p in tables.items():
            check(sorted(p)==list(range(size)),model+' '+name+' is a full permutation')
            hashes[model][name]=hashlib.sha256(','.join(map(str,p)).encode('ascii')).hexdigest()
        if model in EXPECTED_HASHES:check(hashes[model]==EXPECTED_HASHES[model],model+' agrees with frozen previous-stage tables')
        p=[tables['C'][tables['HA'][tables['TA'][i]]] for i in range(size)]
        blocks[model]=p;check(sorted(p)==list(range(size)),model+' mixed block is a permutation')
        cycles[model]=cycle_counts(p);orders[model]=lcm(*cycles[model])
        check(cycles[model]==expected_cycles[model],model+' exact cycle multiplicities')
        check(permutation_power(p,orders[model])==list(range(size)),model+' entire controller returns at its order')
        check(permutation_power(p,96)==list(range(size)),model+' entire controller returns at 96 blocks')
    # Independent exact matrix derivation of the nominal probability formula.
    eye=ident(2);ii=F({16:1});u=root(2)*Q(1,2)
    t=((F(1),F()),(F(),(1+ii)*u));h=scale(((F(1),F(1)),(F(1),F(-1))),u)
    z=((F(1),F()),(F(),F(-1)));x=((F(),F(1)),(F(1),F()))
    cm=tuple(tuple(F(j==((i^1) if i>=2 else i)) for j in range(4)) for i in range(4))
    v=mm(h,t);um=mm(z,v);U=mm(cm,kron(v,eye));xb=kron(eye,x)
    check(mm(U,xb)==mm(xb,U),'mixed unitary conserves target X sectors')
    for sign,small in ((1,v),(-1,um)):
        ket=((u,),(sign*u,));embed=kron(eye,ket)
        check(mm(mm(adj(embed),U),embed)==small,'exact target-X sector restriction')
    vec=tuple(U[i][0] for i in range(4))
    check(vec==(u,F(),F(),u),'one mixed block maps 00 to a Bell state')
    taup=sqabs(trace(v));taum=sqabs(trace(um))
    check(taup==1-u and taum==1+u,'sector squared traces are 1 minus/plus 1/sqrt(2)')
    for tau in (taup,taum):check(2*tau*tau-4*tau+1==0,'sector trace invariant satisfies primitive polynomial 2x^2-4x+1')
    check(taup*taum==Q(1,2) and isqrt(8)**2!=8,'nonintegral field norm and irreducible discriminant for irrational eigenphase proof')
    xp=-(2+root(2))*Q(1,4);xm=-(2-root(2))*Q(1,4)
    cp=[F(1),xp];cmn=[F(1),xm]
    for k in range(2,145):cp.append(2*xp*cp[-1]-cp[-2]);cmn.append(2*xm*cmn[-1]-cmn[-2])
    ideal=[];mat=ident(4)
    for k in range(145):
        prob=((1-cp[k])*(6-root(2))+(1-cmn[k])*(6+root(2)))*Q(1,34)
        check(set(prob.c)<={0},'ideal first-qubit probability is rational by conjugate cancellation')
        ideal.append(prob.c.get(0,Q(0)))
        check(0<=ideal[-1]<=1,'ideal probability lies in the unit interval')
        if k<=96:
            direct=sum((sqabs(mat[i][0]) for i in (2,3)),F())
            check(direct==prob,'full 4x4 matrix power agrees with the Chebyshev formula')
            mat=mm(U,mat)
    check((6-root(2)+6+root(2))*Q(1,34)==Q(6,17),'Cesaro mean coefficient is 6/17')
    check(Q(6,17)-Q(1,3)==Q(1,51),'positive margin for a strict one-third discrepancy')
    check(51*Q(12,17)==36,'explicit universal finite-depth bound coefficient')
    # For small Q, verify the arithmetic used in the analytic root-separation bound.
    for q in range(1,17):
        sp=(1-cp[q])*Q(1,2);sm=(1-cmn[q])*Q(1,2)
        norm=sp*sm
        check(set(norm.c)<={0} and norm.c[0]>0,'small-Q sine-square conjugate norm is positive rational')
        scaled=norm*(4*16**q)
        check(set(scaled.c)=={0} and scaled.c[0].denominator==1 and scaled.c[0]>=1,'integer norm denominator bound (finite examples)')
    def initial(k,l):return enc(k%2,l%2,index[compose(power(SA,k//2),power(SB,l//2))])
    zero=state({3:1,12:1,15:1});rows=[]
    for count in (1,2,4,8,16,24,32,48,72,96,144):
        row={'blocks':count,'elementary_commands':3*count,'ideal_A1':str(ideal[count]),'models':{}}
        for model,p in blocks.items():
            pn=permutation_power(p,count);values=[]
            for k in range(8):
                for l in range(8):
                    a,b,g=dec(pn[initial(k,l)])
                    rr=phase_coeff(act(group[g],phase_coeff(zero,-k,-l)),a,b)
                    prob=(1-rr[12])*Q(1,2)
                    check(set(prob.c)<={0},model+' sampled context probability is rational')
                    value=prob.c.get(0,Q(0));check(0<=value<=1,model+' sampled context probability is valid');values.append(value)
            row['models'][model]={'calibrated_A1':str(values[0]),'uniform_A1':str(sum(values,Q(0))/64),'minimum_context_A1':str(min(values)),'maximum_context_A1':str(max(values))}
            if count%orders[model]==0:check(all(p0==0 for p0 in values),model+' every phase context has zero A1 at a controller return')
        rows.append(row)
    check(ideal[16]==Q(157,256),'original-table return discrepancy at 16 blocks')
    check(ideal[24]==Q(1813,4096),'improved-table return discrepancy at 24 blocks')
    check(ideal[96]==Q(124202613707725,281474976710656)>Q(1,3),'common 96-block discrepancy exceeds one third')
    check(ideal[48]==Q(2421629,16777216)<Q(1,3),'return discrepancies are not uniformly above threshold at every depth')
    TABLE={'schema':'raqm-mixed-depth-v1','chronological_block':['TA','HA','C'],'frames':size,'cycle_counts':cycles,'block_orders':orders,'rows':rows,
      'scope':'Three coarse-table checks. Universal matching independence, all-resolution refinement and Cesaro limits are analytic results.'}
    SUMMARY.update({'frames':size,'cycle_counts':cycles,'block_orders':orders,'common_example_blocks':96,'common_example_commands':288,
      'finite_A1_at_common_example':'0 for every preparation context in all three tables','ideal_A1_at_common_example':str(ideal[96]),
      'universal_return_subsequence_mean':'6/17','uniform_positive_threshold':'1/3',
      'scope':TABLE['scope'],'sha256_gate_permutations':hashes})

EXPECTED_HASHES={'improved': {'HA': 'd28dd9591b8cbf82d8b7782e51d2300257824ec788d868dafcb63a55cdf70332', 'HB': 'b1a43923c6011d37c255968f0f67e87121116ac037d4f9e6d1aaf4be3f5e9298', 'TA': 'f82e962c784534125eb116a9bbc8d5da242607cd2f7c6555e9f3eb84381852f2', 'TB': 'fea55f9c9b1dc459b4e4b0dedb9110226c1f9c6a7d7cd739c740e4018421cc24', 'C': '411bb85521b38374bf50992baf39a641f562b32ae5623937db51999fc44bd0f9'}, 'anchor_preserving': {'HA': '3728a0f3677152b9705b4f2b880c4401b2187f44aae7c092dc6ab1717e1b3a7a', 'HB': '367acd37d052590e206708557ffd5b3ccbe716cb9526056911c5df3b2323dffe', 'TA': 'f82e962c784534125eb116a9bbc8d5da242607cd2f7c6555e9f3eb84381852f2', 'TB': 'fea55f9c9b1dc459b4e4b0dedb9110226c1f9c6a7d7cd739c740e4018421cc24', 'C': 'd972aaba99bfef46dbcbc72d28f75c483154337d838219a9d182c26640d2ae67'}}
if __name__=='__main__':
    started=datetime.now(timezone.utc);out=ROOT/'reproducibility'/('mixed_depth_python_'+started.strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8]);out.mkdir(parents=True)
    status='failed';error=None
    try:main();status='passed'
    except Exception:error=traceback.format_exc();print(error)
    if TABLE:(out/'mixed_depth_results.json').write_text(json.dumps(TABLE,separators=(',',':'))+'\n')
    record={'script':Path(__file__).name,'script_version':VERSION,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      'started_utc':started.isoformat(),'finished_utc':datetime.now(timezone.utc).isoformat(),'python_version':sys.version,'platform':platform.platform(),
      'status':status,'arithmetic':'exact integer and algebraic arithmetic; no floating-point verification','summary':SUMMARY,
      'passed_assertions':sum(c['passed'] for c in CHECKS),'total_assertions':len(CHECKS),'checks':CHECKS,'error':error}
    (out/'mixed_depth_run_record.json').write_text(json.dumps(record,indent=2)+'\n')
    (out/'mixed_depth_assertion_log.txt').write_text('\n'.join(('PASS: ' if c['passed'] else 'FAIL: ')+c['label'] for c in CHECKS)+'\n'+(error or ''))
    print(json.dumps(SUMMARY,indent=2));print('Mathematical verification:',status,'; assertions:',record['passed_assertions']);print('Local record:',out)
    sys.exit(0 if status=='passed' else 1)
