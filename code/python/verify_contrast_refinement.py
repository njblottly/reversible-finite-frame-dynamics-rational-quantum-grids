#!/usr/bin/env python3
"""Contrast robustness and improved matching verifier v0.1. Python 3.10+, standard library.
Run: python3 verify_contrast_refinement.py
Generates the full 46080-frame table and verifies exact probabilities and bounds.
Writes its table, JSON run record and assertion log in a fresh local directory.
Refinement is an analytic theorem in the manuscript; no large refined catalogue is enumerated here.
"""
from collections import deque
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

def main():
    global TABLE
    group,index=enumerate_group(GENS.values());n=len(group)
    check(n==11520,'exact Clifford base size')
    enc=lambda a,b,g:(2*a+b)*n+g
    def dec(x):q,g=divmod(x,n);return q//2,q%2,g
    left={name:[index[compose(p,g)] for g in group] for name,p in GENS.items()}
    twist=compose(inv(SB),compose(C,SB))
    twisted=[index[compose(twist,g)] for g in group]
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
    hashes={}
    for model,tables in [('improved',new),('anchor_preserving',anchored)]:
        hashes[model]={}
        for name,p in tables.items():
            check(sorted(p)==list(range(4*n)),model+' '+name+' full permutation')
            hashes[model][name]=hashlib.sha256(','.join(map(str,p)).encode('ascii')).hexdigest()
            if name in ('HA','HB','C'):check(all(p[p[i]]==i for i in range(4*n)),model+' '+name+' involution')
            else:check(all(p[p[p[p[p[p[p[p[i]]]]]]]]==i for i in range(4*n)),model+' '+name+' eighth power identity')
        for ga,gb in [('HA','HB'),('HA','TB'),('TA','HB'),('TA','TB')]:
            check(all(tables[ga][tables[gb][i]]==tables[gb][tables[ga][i]] for i in range(4*n)),model+' opposite local commands commute')
    if EXPECTED_HASHES:check(hashes==EXPECTED_HASHES,'all frozen permutation hashes agree')
    for g in (ID,C,HA,compose(C,HA)):
        x=enc(0,0,index[g]);check(anchored['C'][x]==enc(0,0,index[compose(C,g)]),'anchor-preserving C table retains exact calibration')
    # Exact algebraic certificates: all relative spectra consist of 1 and zeta or its conjugate.
    eye=ident(2);eye4=ident(4);ii=F({16:1});u=root(2)*Q(1,2);zeta=(1+ii)*u
    t=((F(1),F()),(F(),zeta));h=scale(((F(1),F(1)),(F(1),F(-1))),u);s=mm(t,t)
    cm=tuple(tuple(F(j==((i^1) if i>=2 else i)) for j in range(4)) for i in range(4))
    h4=kron(h,eye)
    for a in (0,1):
        actual=mm(t,h) if a==0 else mm(h,adj(t));rel=mm(adj(h),actual);z=zeta if a==0 else zeta.conj()
        check(mm(addmat(rel,scale(eye,-1)),addmat(rel,scale(eye,-z)))==scale(eye,0),'improved H relative spectral polynomial')
        check(trace(rel)==1+z,'improved H has both endpoint eigenvalues')
    for a in (0,1):
        for b in (0,1):
            before=kron(t if a else eye,t if b else eye);after=kron(t if a else eye,t if not b else eye)
            actual=mm(mm(after,cm),adj(before));rel=mm(adj(cm),actual);z=zeta if b==0 else zeta.conj()
            check(mm(addmat(rel,scale(eye4,-1)),addmat(rel,scale(eye4,-z)))==scale(eye4,0),'improved CNOT relative spectral polynomial')
            check(trace(rel)==2+2*z,'improved CNOT spectrum has multiplicities two and two')
            # Anchored twist gives the opposite local conjugation for target bit one.
            e=before;mid=cm if b==0 else mm(mm(kron(eye,adj(s)),cm),kron(eye,s))
            actual=mm(mm(e,mid),adj(e));rel=mm(adj(cm),actual);dd=addmat(eye4,scale(rel,-1))
            expected=tuple(tuple((2-root(2) if b and i==j and i>=2 else F()) for j in range(4)) for i in range(4))
            check(mm(adj(dd),dd)==expected,'anchor-preserving C error unchanged and certified')
    check(2+root(2)-(1+u)*(1+u)==Q(1,2),'improved error strictly below previous H and C bounds')
    def initial(k,l):return enc(k%2,l%2,index[compose(power(SA,k//2),power(SB,l//2))])
    def end(x,word,tables):
        for gate in word:x=tables[gate][x]
        return x
    def direct(word,k,l,r,tables):
        a,b,g=dec(end(initial(k,l),word,tables))
        return phase_coeff(act(group[g],phase_coeff(r,-k,-l)),a,b)
    zero=state({3:1,12:1,15:1});phi=state({5:1,10:-1,15:1});pp=state({1:1,4:1,5:1})
    rows=[]
    for power_n in range(9):
        word=['HA','C']+['TB']*power_n+['C','HA'];values=[]
        for k in range(8):
            for l in range(8):
                p=computational(direct(word,k,l,zero,new))[0]
                pa=computational(direct(word,k,l,zero,anchored))[0]
                old_flipped=computational(direct(word,k,l^1,zero,old))[0]
                check(p==pa==old_flipped,'both new fringes exchange the original target-parity branches')
                values.append(p)
        avg=sum(values,F())*Q(1,64)
        rows.append({'n':power_n,'calibrated_p00':encode_scalar(values[0]),'uniform_p00':encode_scalar(avg)})
    for model,tables in [('improved',new),('anchor_preserving',anchored)]:
        for k in range(8):
            for l in range(8):
                p1=computational(direct(['HA','C','TB','C','HA'],k,l,zero,tables))[0]
                p3=computational(direct(['HA','C','TB','TB','TB','C','HA'],k,l,zero,tables))[0]
                check(p1-p3==Q(1,2),model+' contrast remains 1/2 in every context')
        # Common channel is verified independently through the cancellation of even initial phase parts.
        for word in (['TA','HA','C','TB','HB','C'],['HB','C','TA','HA','TB','C','HB']):
            for k in range(8):
                for l in range(8):
                    for r in (zero,phi,pp):check(direct(word,k,l,r,tables)==direct(word,k%2,l%2,r,tables),model+' preparation depends only on phase parities')
    # Exact rational comparison for a sufficient refinement threshold, not measured refined-table errors.
    check(Q(1,2)>Q(5,8)**2,'if both gate errors <=1/64, contrast lower bound exceeds 1/2')
    TABLE={'schema':'raqm-contrast-matching-variants-v1','clifford_signed_pauli_actions':group,
      'label_encoding':'zero-based (2*a+b)*11520+g_index','improved_permutations':new,
      'anchor_preserving_C_permutation':anchor_c,'canonical_permutation_hashes':hashes,
      'note':'Anchor-preserving variant uses original HA/HB and the unchanged exact phase tables. Genuine refinement is proved analytically, not enumerated here.'}
    SUMMARY.update({'frames':4*n,'improved_forward_edges':5*4*n,'hashes':hashes,
      'improved_sharp_H_and_C_error':'2*sin(pi/16)=sqrt(2-sqrt(2+sqrt(2)))',
      'anchor_variant_H_error':'sqrt(1-1/sqrt(2))','anchor_variant_C_error':'sqrt(2-sqrt(2))',
      'fringe_rows':rows,'contrast_all_tested_models':'1/2',
      'refinement_theorem':'For errors epsH,epsC and exact phase steps: abs(contrast-1/sqrt(2)) <= 4*(epsH+epsC), uniformly over fixed priors.',
      'scope':'Finite checks certify same-catalogue variants. Nested catalogue refinement and uniform convergence are analytic results, not numerically instantiated catalogues.'})

EXPECTED_HASHES={'improved': {'HA': 'd28dd9591b8cbf82d8b7782e51d2300257824ec788d868dafcb63a55cdf70332', 'HB': 'b1a43923c6011d37c255968f0f67e87121116ac037d4f9e6d1aaf4be3f5e9298', 'TA': 'f82e962c784534125eb116a9bbc8d5da242607cd2f7c6555e9f3eb84381852f2', 'TB': 'fea55f9c9b1dc459b4e4b0dedb9110226c1f9c6a7d7cd739c740e4018421cc24', 'C': '411bb85521b38374bf50992baf39a641f562b32ae5623937db51999fc44bd0f9'}, 'anchor_preserving': {'HA': '3728a0f3677152b9705b4f2b880c4401b2187f44aae7c092dc6ab1717e1b3a7a', 'HB': '367acd37d052590e206708557ffd5b3ccbe716cb9526056911c5df3b2323dffe', 'TA': 'f82e962c784534125eb116a9bbc8d5da242607cd2f7c6555e9f3eb84381852f2', 'TB': 'fea55f9c9b1dc459b4e4b0dedb9110226c1f9c6a7d7cd739c740e4018421cc24', 'C': 'd972aaba99bfef46dbcbc72d28f75c483154337d838219a9d182c26640d2ae67'}}
if __name__=='__main__':
    started=datetime.now(timezone.utc);out=ROOT/'reproducibility'/('contrast_python_'+started.strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8]);out.mkdir(parents=True)
    status='failed';error=None
    try:main();status='passed'
    except Exception:error=traceback.format_exc();print(error)
    if TABLE:(out/'contrast_matching_tables.json').write_text(json.dumps(TABLE,separators=(',',':'))+'\n')
    record={'script':Path(__file__).name,'script_version':VERSION,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      'started_utc':started.isoformat(),'finished_utc':datetime.now(timezone.utc).isoformat(),'python_version':sys.version,'platform':platform.platform(),
      'status':status,'arithmetic':'exact integer and algebraic arithmetic; no floating-point verification','summary':SUMMARY,
      'passed_assertions':sum(c['passed'] for c in CHECKS),'total_assertions':len(CHECKS),'checks':CHECKS,'error':error}
    (out/'contrast_refinement_run_record.json').write_text(json.dumps(record,indent=2)+'\n')
    (out/'contrast_refinement_assertion_log.txt').write_text('\n'.join(('PASS: ' if c['passed'] else 'FAIL: ')+c['label'] for c in CHECKS)+'\n'+(error or ''))
    print(json.dumps(SUMMARY,indent=2));print('Mathematical verification:',status,'; assertions:',record['passed_assertions']);print('Local record:',out)
    sys.exit(0 if status=='passed' else 1)
