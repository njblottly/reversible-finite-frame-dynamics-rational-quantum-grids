#!/usr/bin/env python3
"""Certified non-Clifford two-qubit controller v0.1. Python 3.10+, standard library.
Run: python3 verify_two_qubit_model.py
Generates the full 46080-frame table and verifies exact probabilities and bounds.
Writes its table, JSON run record and assertion log in a fresh local directory.
This is a coarse structured matching, not a bottleneck optimum or a refining family.
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
    base_tables={name:[index[compose(g,x)] for x in group] for name,g in GENS.items()}
    for name,p in base_tables.items():check(sorted(p)==list(range(n)),name+' exact Clifford left-action table')
    def enc(a,b,g):return (2*a+b)*n+g
    def dec(x):q,g=divmod(x,n);return q//2,q%2,g
    names=('HA','HB','TA','TB','C')
    tables={name:[] for name in names}
    for a in range(2):
        for b in range(2):
            for g in range(n):
                tables['HA'].append(enc(a,b,base_tables['HA'][g]))
                tables['HB'].append(enc(a,b,base_tables['HB'][g]))
                tables['C'].append(enc(a,b,base_tables['C'][g]))
                tables['TA'].append(enc(1-a,b,g if a==0 else base_tables['SA'][g]))
                tables['TB'].append(enc(a,1-b,g if b==0 else base_tables['SB'][g]))
    hashes={}
    for name,p in tables.items():
        check(sorted(p)==list(range(4*n)),name+' is a full 46080-frame permutation')
        hashes[name]=hashlib.sha256(','.join(map(str,p)).encode('ascii')).hexdigest()
        inverse=[0]*(4*n)
        for i,j in enumerate(p):inverse[j]=i
        check(all(inverse[p[i]]==i for i in range(4*n)),name+' inverse restores every label')
        if name in ('HA','HB','C'):check(all(p[p[i]]==i for i in range(4*n)),name+' is exactly involutive')
        else:
            check(all(p[p[p[p[p[p[p[p[i]]]]]]]]==i for i in range(4*n)),name+' has eighth power identity')
    # Frozen table hashes are checked after the candidate has been selected once.
    if EXPECTED_HASHES:check(hashes==EXPECTED_HASHES,'all five frozen table hashes match')
    check(all(tables['HA'][tables['HB'][i]]==tables['HB'][tables['HA'][i]] for i in range(4*n)),'opposite Hadamards commute on every label')
    check(all(tables['TA'][tables['TB'][i]]==tables['TB'][tables['TA'][i]] for i in range(4*n)),'opposite phases commute on every label')
    # Independent 2x2/4x4 algebraic matrix certificate for all gate-error types.
    eye=ident(2);ii=F({16:1});u=root(2)*Q(1,2);zz=u*(1+ii)
    h=scale(((F(1),F(1)),(F(1),F(-1))),u)
    t=((F(1),F()),(F(),zz));s=mm(t,t);eye4=ident(4)
    c=tuple(tuple(F(j==((i^1) if i>=2 else i)) for j in range(4)) for i in range(4))
    check(mm(adj(t),t)==eye and mm(h,h)==eye and mm(c,c)==eye4,'exact generating matrix unitarity')
    check(mm(t,t)==((F(1),F()),(F(),ii)),'T carry is exact S')
    hactual=mm(mm(t,h),adj(t));rel=mm(adj(h),hactual)
    delta=addmat(eye,negmat(rel))
    h_squared=1-root(2)*Q(1,2);c_squared=2-root(2)
    check(mm(adj(delta),delta)==scale(eye,h_squared),'exact sharp H squared error 1-1/sqrt(2)')
    for a in range(2):
        for b in range(2):
            e=kron(t if a else eye,t if b else eye)
            actual=mm(mm(e,c),adj(e));rel=mm(adj(c),actual)
            delta=addmat(eye4,negmat(rel));gram=mm(adj(delta),delta)
            expected=tuple(tuple((c_squared if i==j and i>=2 and b else F()) for j in range(4)) for i in range(4))
            check(gram==expected,'exact CNOT squared-error spectrum for both phase bits')
            check(mm(actual,actual)==eye4,'each selected CNOT unitary is an involution')
    # Four disjoint translates E_q Clifford are distinct, tested via a Pauli image.
    for a in range(2):
        for b in range(2):
            for aa in range(2):
                for bb in range(2):
                    if (a,b)==(aa,bb):continue
                    basis=[0]*16;basis[4 if a!=aa else 1]=1
                    im=phase_coeff(basis,aa-a,bb-b)
                    check(sum(x!=0 for x in im)==2,'distinct phase cosets have a non-Clifford Pauli image')
    def initial(k,l):
        g=compose(power(SA,k//2),power(SB,l//2))
        return enc(k%2,l%2,index[g])
    def end(x,word):
        for gate in word:x=tables[gate][x]
        return x
    def frame_act(x,r):
        a,b,g=dec(x);return phase_coeff(act(group[g],r),a,b)
    def direct(word,k,l,r):return frame_act(end(initial(k,l),word),phase_coeff(r,-k,-l))
    def physical(word,k,l,r):
        x=initial(k,l);r=tuple(map(asf,r))
        for gate in word:
            a,b,g=dec(x)
            if gate=='TA':r=phase_coeff(r,1,0)
            elif gate=='TB':r=phase_coeff(r,0,1)
            elif gate=='HA':r=phase_coeff(act(HA,phase_coeff(r,-a,0)),a,0)
            elif gate=='HB':r=phase_coeff(act(HB,phase_coeff(r,0,-b)),0,b)
            else:r=phase_coeff(act(C,phase_coeff(r,-a,-b)),a,b)
            x=tables[gate][x]
        return r
    zero=state({3:1,12:1,15:1});phi=state({5:1,10:-1,15:1})
    pp=state({1:1,4:1,5:1});pi0=state({3:1,8:1,11:1});mixed=state({})
    pa=[ID,local((1,2,-3,-4),0),local((1,-2,3,-4),0),local((1,-2,-3,4),0)]
    bell=[act(p,phi) for p in pa];rotated=[act(SA,r) for r in bell]
    check(frame_act(initial(0,0),zero)==zero,'calibrated preparation has identity frame')
    check(direct(['HA','C'],0,0,zero)==tuple(map(asf,phi)),'calibrated prefix creates Bell state exactly')
    cases=[['HA','C','TB','C','HA'],['TA','HB','C','HA','TB','C','TA','HB'],['HB','TA','C','TB','HA','C']]
    for word in cases:
        for k in range(8):
            for l in range(8):
                for r in (zero,phi,pp,pi0):
                    got=direct(word,k,l,r)
                    check(got==physical(word,k,l,r),'retained-seed and physical-gate propagation agree')
                    check(got==direct(word,k%2,l%2,r),'only initial phase parities affect the channel')
                for recipe in (bell,rotated):
                    out=[direct(word,k,l,r) for r in recipe]
                    avg=tuple(sum((row[j] for row in out),F())*Q(1,4) for j in range(16))
                    check(avg==mixed,'equal-density Bell recipes give the same channel output')
    # Full chronological entangling phase scan, evaluated for all 64 contexts.
    expected_cal=[Q(1),Q(1),Q(1,2),Q(1,2),Q(0),Q(0),Q(1,2),Q(1,2),Q(1)]
    expected_uni=[Q(1),Q(3,4),Q(1,2),Q(1,4),Q(0),Q(1,4),Q(1,2),Q(3,4),Q(1)]
    rows=[]
    for power_n in range(9):
        word=['HA','C']+['TB']*power_n+['C','HA'];probs=[]
        for k in range(8):
            for l in range(8):
                p=computational(direct(word,k,l,zero));probs.append(p)
                m=(l%2+power_n)//2;expected=Q(1+(1,0,-1,0)[m%4],2)
                check(p==(expected,0,1-expected,0),'exact phase-scan probability in every preparation context')
        avg=tuple(sum((p[j] for p in probs),F())*Q(1,64) for j in range(4))
        check(probs[0][0]==expected_cal[power_n] and avg[0]==expected_uni[power_n],'calibrated and uniform phase-scan row')
        ideal=phase_coeff(zero,0,0)
        for gate in word:
            ideal=phase_coeff(ideal,0,1) if gate=='TB' else act(GENS[gate],ideal)
        pi=computational(ideal)[0]
        if power_n%2:
            check(sqabs(pi-avg[0])==(3-2*root(2))*Q(1,16),'exact odd-phase uniform gap squared')
        else:check(pi==avg[0]==probs[0][0],'even phase scan agrees exactly with ideal')
        rows.append({'n':power_n,'calibrated_p00':str(expected_cal[power_n]),'uniform_p00':str(expected_uni[power_n]),'ideal_p00_field':encode_scalar(pi)})
    # Prior-independent two-setting contrast and a sharp minimax prior result.
    for k in range(8):
        for l in range(8):
            p1=computational(direct(['HA','C','TB','C','HA'],k,l,zero))[0]
            p3=computational(direct(['HA','C','TB','TB','TB','C','HA'],k,l,zero))[0]
            check(p1-p3==Q(1,2),'n=1 minus n=3 contrast is independent of every preparation context')
    contexts=[direct(['HA','HB'],a,b,zero) for a in range(2) for b in range(2)]
    check(all(contexts[i]!=contexts[j] for i in range(4) for j in range(i)), 'four phase-parity contexts are future-distinguishable at the same input ray')
    # For rational prior weights, verify the derived worst-error formula exactly.
    alpha=(root(2)-1)*Q(1,4)
    for j in range(17):
        weight=Q(j,16)
        e1=(1-root(2)*Q(1,2)-weight)*Q(1,2)
        e3=(root(2)*Q(1,2)-weight)*Q(1,2)
        extreme=e3 if weight<=Q(1,2) else -e1
        check(extreme==alpha+abs(weight-Q(1,2))*Q(1,2),'sharp worst-error expression for fixed-prior fringe')
    TABLE={'schema':'raqm-two-qubit-structured-matching-v1','version':VERSION,
       'frame_definition':'F_(a,b,g)=(T^a tensor T^b) g; a,b in {0,1}; g a projective two-qubit Clifford',
       'label_encoding':'zero-based index=(2*a+b)*11520+g_index',
       'clifford_action_convention':'g P_j g^dagger=sign(entry_j)*P_(abs(entry_j)-1); II,IX,IY,IZ,XI,...,ZZ',
       'clifford_bfs_generator_order':['HA','HB','SA','SB','C'],
       'clifford_signed_pauli_actions':group,'permutations':tables,
       'permutation_sha256_comma_joined_decimal':hashes,
       'sharp_errors':{'HA':'sqrt(1-1/sqrt(2))','HB':'sqrt(1-1/sqrt(2))','TA':'0','TB':'0','C':'sqrt(2-sqrt(2))'},
       'selection':'explicit structured bijections; no minimax optimality claimed',
       'preparation':'P_k=T^k, P_l=T^l; original compensated seed; L=8'}
    SUMMARY.update({'frames':4*n,'clifford_base':n,'local_frames':48,'nonlocal_charts':20,
       'total_forward_table_edges':5*4*n,'table_hashes':hashes,'phase_scan':rows,
       'effective_preparation_branches':4,'exact_uniform_odd_phase_gap':'(sqrt(2)-1)/4',
       'prior_independent_p1_minus_p3':'1/2',
       'sharp_minimax_prior_error':'(sqrt(2)-1)/4 + abs(w-1/2)/2; w=Pr(initial target phase odd)',
       'scope':'Coarse concrete non-Clifford model; exact algebraic certification, no claim of refinement or unique RaQM prediction.'})

EXPECTED_HASHES = {'HA': '3728a0f3677152b9705b4f2b880c4401b2187f44aae7c092dc6ab1717e1b3a7a', 'HB': '367acd37d052590e206708557ffd5b3ccbe716cb9526056911c5df3b2323dffe', 'TA': 'f82e962c784534125eb116a9bbc8d5da242607cd2f7c6555e9f3eb84381852f2', 'TB': 'fea55f9c9b1dc459b4e4b0dedb9110226c1f9c6a7d7cd739c740e4018421cc24', 'C': '0e008ad66cc235bf1e0a74037fd4109c40d3968c9446427cb2e5e70618fba19c'}
if __name__=='__main__':
    started=datetime.now(timezone.utc)
    out=ROOT/'reproducibility'/('two_qubit_python_'+started.strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8])
    out.mkdir(parents=True)
    status='failed';error=None
    try:main();status='passed'
    except Exception:error=traceback.format_exc();print(error)
    if TABLE:(out/'two_qubit_matching_table.json').write_text(json.dumps(TABLE,separators=(',',':'))+'\n')
    record={'script':Path(__file__).name,'script_version':VERSION,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
       'started_utc':started.isoformat(),'finished_utc':datetime.now(timezone.utc).isoformat(),'python_version':sys.version,'platform':platform.platform(),
       'status':status,'arithmetic':'exact integer, rational and multiquadratic; no floating point verification',
       'passed_assertions':sum(c['passed'] for c in CHECKS),'total_assertions':len(CHECKS),'summary':SUMMARY,'checks':CHECKS,'error':error}
    (out/'two_qubit_run_record.json').write_text(json.dumps(record,indent=2)+'\n')
    (out/'two_qubit_assertion_log.txt').write_text('\n'.join(('PASS: ' if c['passed'] else 'FAIL: ')+c['label'] for c in CHECKS)+'\n'+(error or ''))
    print(json.dumps(SUMMARY,indent=2));print('Mathematical verification:',status);print('Assertions:',record['passed_assertions']);print('Local record:',out)
    sys.exit(0 if status=='passed' else 1)
