#!/usr/bin/env python3
"""Two-qubit phase lift and whole-history certificates. Standalone Python 3.10+.
Standard library only. Exact signed Pauli actions, algebraic matrices and
rational probability comparisons. Records and hashes saved beside this file.
The 184320 objects are retained controller labels, not distinct matrices.
The seven-block obstruction covers the classified local Clifford-covariant
Hadamard family with the specified CNOT table, not all allowed matchings.
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


def bounds(x,bits=100):
    lo=hi=Q(0);m=2**bits
    assert all(not k&16 for k in x.c)
    for k,c in x.c.items():
        z=FACTORS[k];r=isqrt(z*m*m);a=Q(r,m);b=a if r*r==z*m*m else Q(r+1,m)
        lo+=c*(a if c>=0 else b);hi+=c*(b if c>=0 else a)
    return lo,hi

def sign(x):
    if x==0:return 0
    for bits in (100,200,400,800):
        lo,hi=bounds(x,bits)
        if lo>0:return 1
        if hi<0:return -1
    raise ArithmeticError('Exact sign unresolved')

def main():
    global TABLE
    eye=ident(2);ii=F({16:1});u=root(2)*Q(1,2)
    h=scale(((F(1),F(1)),(F(1),F(-1))),u)
    t=((F(1),F()),(F(),(1+ii)*u));s=mm(t,t)
    x=((F(),F(1)),(F(1),F()));y=((F(),-ii),(ii,F()));z=((F(1),F()),(F(),F(-1)))
    cm=tuple(tuple(F(j==((i^1) if i>=2 else i)) for j in range(4)) for i in range(4))
    paulis=[kron(a,b) for a in (eye,x,y,z) for b in (eye,x,y,z)]
    # Exact nonclosure witness C T_B. All four possible left phase translates
    # would have to leave a Clifford factor, but its X_B image is not Pauli.
    supports=[]
    for a in range(2):
        for b in range(2):
            E=kron(t if a else eye,t if b else eye)
            K=mm(mm(adj(E),cm),kron(eye,t))
            image=mm(mm(K,paulis[1]),adj(K))
            coeff=[trace(mm(P,image))*Q(1,4) for P in paulis]
            support=[j for j,c in enumerate(coeff) if c!=0]
            check(len(support)==(2 if b==0 else 4),'C T_B absent from each original phase-Clifford translate')
            supports.append({'a':a,'b':b,'Pauli_support':support,'coefficients':[encode_scalar(coeff[j]) for j in support]})
    for a in range(2):
        for b in range(2):
            E=kron(t if a else eye,t if b else eye);rel=mm(adj(cm),mm(mm(E,cm),adj(E)))
            diag=(F(1),F(1),(1+ii)*u,(1-ii)*u) if b else (F(1),)*4
            check(rel==tuple(tuple(diag[i] if i==j else F() for j in range(4)) for i in range(4)), 'exact CNOT relative spectrum certifies its unchanged error bound')
    # Complete one-qubit Clifford classification in exact matrix arithmetic.
    lid=(1,2,3,4);local_actions=[lid];local_data={lid:(eye,'')}
    for g in local_actions:
        for name,action,mat in [('H',h1,h),('S',s1,s)]:
            k=compose(action,g)
            if k not in local_data:
                local_data[k]=(mm(mat,local_data[g][0]),local_data[g][1]+name);local_actions.append(k)
    check(len(local_actions)==24,'complete local Clifford group')
    threshold=F(Q(3,2))+root(2);candidates=[]
    for g in local_actions:
        mat,word=local_data[g];physical=mm(mm(t,mat),adj(t))
        q=sqabs(trace(mm(adj(h),physical)))
        if sign(q-threshold)>=0:candidates.append((g,mat,word))
    check([v[2] for v in candidates]==['H','SH','HSSS','SHSSS'],'exactly four odd-parity Hadamard choices meet the original bound')
    expected_matrices=[h,mm(h,s),mm(mm(mm(adj(s),eye),h),eye),mm(mm(adj(s),h),s)]
    # Word strings are chronological; SH means H S, HSSS means S^-1 H.
    check(all(sqabs(trace(mm(adj(a[1]),b)))==4 for a,b in zip(candidates,expected_matrices)),'closed-form matrices for the four choices')
    for g,mat,word in candidates:
        check(sqabs(trace(mm(adj(h),mm(mm(t,mat),adj(t)))))==threshold,'all four choices attain the sharp Hadamard error bound')
    check(sqabs(trace(mm(adj(mm(t,h)),mm(h,t))))==threshold,'duplicate physical frame labels have distinct H successors in the inherited table')
    group,index=enumerate_group(GENS.values());n=len(group);M=16*n
    check(n==11520 and M==184320,'lift has 184320 retained labels')
    actions={'HA0':HA,'HB0':HB,'SA':SA,'SB':SB,'C':C}
    for k,(g,mat,word) in enumerate(candidates):actions['HA'+str(k+1)]=local(g,0);actions['HB'+str(k+1)]=local(g,1)
    left={name:[index[compose(p,g)] for g in group] for name,p in actions.items()}
    rightA=[index[compose(g,SA)] for g in group];rightB=[index[compose(g,SB)] for g in group]
    enc=lambda a,b,g,r,s:(((2*a+b)*n+g)*4+2*r+s)
    def dec(i):q,rs=divmod(i,4);ab,g=divmod(q,n);return ab//2,ab%2,g,rs//2,rs%2
    from array import array
    rhoA=array('I');rhoB=array('I');TA=array('I');TB=array('I');PC=array('I')
    HAs=[array('I') for _ in candidates];HBs=[array('I') for _ in candidates]
    for a in range(2):
        for b in range(2):
            for g in range(n):
                for r in range(2):
                    for ss in range(2):
                        rhoA.append(enc(a,b,rightA[g] if r else g,1-r,ss))
                        rhoB.append(enc(a,b,rightB[g] if ss else g,r,1-ss))
                        TA.append(enc(1-a,b,left['SA'][g] if a else g,r,ss))
                        TB.append(enc(a,1-b,left['SB'][g] if b else g,r,ss))
                        PC.append(enc(a,b,left['C'][g],r,ss))
                        for k in range(4):
                            HAs[k].append(enc(a,b,left['HA'+str(k+1) if a else 'HA0'][g],r,ss))
                            HBs[k].append(enc(a,b,left['HB'+str(k+1) if b else 'HB0'][g],r,ss))
    tables={'TA':TA,'TB':TB,'C':PC,**{'HA_'+candidates[k][2]:HAs[k] for k in range(4)},**{'HB_'+candidates[k][2]:HBs[k] for k in range(4)}}
    gate_hashes={}
    for name,p in tables.items():
        check(sorted(p)==list(range(M)),name+' full lifted table is a permutation')
        check(all(p[rhoA[i]]==rhoA[p[i]] and p[rhoB[i]]==rhoB[p[i]] for i in range(M)),name+' commutes with both right-phase actions on every retained label')
        gate_hashes[name]=hashlib.sha256(','.join(map(str,p)).encode()).hexdigest()
    check(all(rhoA[rhoB[i]]==rhoB[rhoA[i]] for i in range(M)),'two right-phase generators commute everywhere')
    for p in (rhoA,rhoB):
        q=list(range(M))
        for _ in range(8):q=[p[i] for i in q]
        check(q==list(range(M)),'right-phase generator has eighth power identity')
    for pa in [TA]+HAs:
        for pb in [TB]+HBs:check(all(pa[pb[i]]==pb[pa[i]] for i in range(M)),'every A-local and B-local table commutes on all labels')
    # Orbit counting; group action is free and all quotient labels retained.
    seen=bytearray(M);orbits=0
    for i in range(M):
        if seen[i]:continue
        orbitset=set();ia=i
        for _ in range(8):
            ib=ia
            for _ in range(8):orbitset.add(ib);ib=rhoB[ib]
            ia=rhoA[ia]
        check(len(orbitset)==64,'right-phase orbit has exactly 64 labels')
        for j in orbitset:seen[j]=1
        orbits+=1
    check(orbits==2880,'right-phase quotient has 2880 labels')
    origin=enc(0,0,index[ID],0,0);contexts=[];ia=origin
    for k in range(8):
        ib=ia
        for l in range(8):contexts.append(ib);ib=rhoB[ib]
        ia=rhoA[ia]
    for p in HAs+HBs:
        check(p[p[origin]]==origin,'Hadamard calibration round trip at identity')
    check(PC[PC[origin]]==origin,'CNOT calibration round trip at identity')
    check(all(PC[PC[p[origin]]]==p[origin] for p in HAs),'CNOT calibrated at H_A and C H_A')
    # Exact ideal mixed-control probabilities through 64 blocks.
    xp=-(2+root(2))*Q(1,4);xm=-(2-root(2))*Q(1,4);cp=[F(1),xp];cn=[F(1),xm]
    for k in range(2,65):cp.append(2*xp*cp[-1]-cp[-2]);cn.append(2*xm*cn[-1]-cn[-2])
    ideal=[];U=mm(cm,kron(mm(h,t),eye));mat=ident(4)
    for k in range(65):
        prob=((1-cp[k])*(6-root(2))+(1-cn[k])*(6+root(2)))*Q(1,34)
        check(set(prob.c)<={0},'ideal entangling probability is exact rational')
        check(prob==sum((sqabs(mat[i][0]) for i in (2,3)),F()),'direct 4x4 ideal power matches formula')
        ideal.append(prob.c.get(0,Q(0)));mat=mm(U,mat)
    zero=state({3:1,12:1,15:1});rows=[]
    for ci,(ga,ha,word) in enumerate(candidates):
        pH=HAs[ci];ends=contexts[:];actual=[];distribution=[];first=None
        for depth in range(65):
            probs=[];dists=[]
            for j,label in enumerate(ends):
                a,b,g,r,ss=dec(label);k,l=divmod(j,8)
                seed=phase_coeff(zero,-k,-l);inner=phase_coeff(seed,r,ss)
                rr=phase_coeff(act(group[g],inner),a,b)
                prob=(1-rr[12])*Q(1,2);check(set(prob.c)<={0},'lifted entangling probability is rational')
                probs.append(prob.c.get(0,Q(0)));dists.append(computational(rr))
            check(all(x==probs[0] for x in probs) and all(x==dists[0] for x in dists),'all 64 contexts have identical joint probabilities at each depth')
            actual.append(probs[0]);distribution.append([encode_scalar(x) for x in dists[0]])
            if first is None and abs(probs[0]-ideal[depth])>Q(1,3):first=depth
            ends=[PC[pH[TA[i]]] for i in ends]
        check(first==[7,6,7,7][ci],'exact first discrepancy greater than one third for '+word)
        check(actual[7]==Q(1,2) and ideal[7]==Q(1,16),'common seven-block gap is 7/16')
        rows.append({'chronological_odd_H_word':word,'matrix_form':['H','H S','S^-1 H','S^-1 H S'][ci],'first_gap_above_1_over_3':first,'A1_probabilities':[str(x) for x in actual],'joint_probabilities':distribution})
    check(all(abs(Q(rows[0]['A1_probabilities'][k])-ideal[k])<=Q(1,3) for k in range(7)),'one fixed table survives every depth through six')
    # Bell preparation is exact for every context when command order is H_A,C.
    end=PC[HAs[0][origin]];a,b,g,r,ss=dec(end);rr=phase_coeff(act(group[g],phase_coeff(zero,r,ss)),a,b)
    check(rr==state({5:1,10:-1,15:1}),'calibrated H_A,C prepares an exact Bell state')
    TABLE={'scope':'Two-bit labelled phase lift of the 46080-frame model; exact 16-member independent local-H family (four distinct repeated-block histories), original C table fixed; no unrestricted matching or refinement optimality claim.',
      'base_frames':46080,'retained_labels':M,'phase_orbit_size':64,'phase_quotient_labels':orbits,'additional_label_bits':2,'labels_need_not_be_distinct_matrices':True,
      'nonclosure_witness':supports,'classified_H1_choices':[r['matrix_form'] for r in rows],
      'H_error_bound':'sqrt(1-1/sqrt(2))','C_error_bound':'sqrt(2-sqrt(2))','ideal_A1_probabilities':[str(x) for x in ideal],'history_rows':rows,
      'common_point_depth_above_one_third':7,'D_curve':8,'common_point_gap':'7/16','gate_table_sha256':gate_hashes}
    SUMMARY.update({k:TABLE[k] for k in ('base_frames','retained_labels','phase_quotient_labels','additional_label_bits','classified_H1_choices','common_point_depth_above_one_third','D_curve','common_point_gap','scope')})

if __name__=='__main__':
    started=datetime.now(timezone.utc);out=ROOT/'reproducibility'/('twoqubit_histories_python_'+started.strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8]);out.mkdir(parents=True)
    status='failed';error=None
    try:main();status='passed'
    except Exception:error=traceback.format_exc();print(error)
    if TABLE:
        (out/'twoqubit_history_results.json').write_text(json.dumps(TABLE,indent=2)+'\n')
        (ROOT/'twoqubit_history_results.json').write_text(json.dumps(TABLE,indent=2)+'\n')
    record={'script':Path(__file__).name,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'started_utc':started.isoformat(),'finished_utc':datetime.now(timezone.utc).isoformat(),'python_version':sys.version,'platform':platform.platform(),'status':status,'arithmetic':'exact algebraic, rational and integer; no floating-point certification','summary':SUMMARY,'total_assertions':len(CHECKS),'passed_assertions':sum(c['passed'] for c in CHECKS),'checks':CHECKS,'error':error}
    (out/'twoqubit_histories_run_record.json').write_text(json.dumps(record,indent=2)+'\n')
    (out/'twoqubit_histories_assertion_log.txt').write_text('\n'.join(('PASS: ' if c['passed'] else 'FAIL: ')+c['label'] for c in CHECKS)+'\n'+(error or ''))
    print(json.dumps(SUMMARY,indent=2));print('Mathematical verification:',status,'Assertions:',len(CHECKS));print('Local record:',out)
    sys.exit(0 if status=='passed' else 1)
