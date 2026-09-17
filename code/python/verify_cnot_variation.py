#!/usr/bin/env python3
"""CNOT matching counterexample: standalone Python 3.10+, standard library.
Run this saved file in PyCharm or: python3 verify_cnot_variation.py
No input files, third-party packages or large matching search are required.
Fresh JSON and text records are written to reproducibility/ beside this file.
This script has been statically checked by the authoring assistant but NOT
executed there. Mathematical runtime verification is left to the user's Mac.
The complete table is certified algebraically by its phase-type update rules;
184320 individual label entries are not enumerated.
"""
from fractions import Fraction as Q
from pathlib import Path
from datetime import datetime, timezone
from math import isqrt
import hashlib,json,platform,traceback,uuid
SCAN_HORIZON=64  # Increase if desired; seven blocks suffice for the counterexample.
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
    source=Path(__file__).resolve();now=lambda:datetime.now(timezone.utc).isoformat()
    folder=source.parent/'reproducibility'/('cnot_variation_python_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8]);folder.mkdir(parents=True)
    started=now();checks=[];results={};status='passed';error=None
    def check(ok,label):
        checks.append({'label':label,'passed':bool(ok)})
        if not ok:raise AssertionError(label)
    try:
        check(isinstance(SCAN_HORIZON,int) and SCAN_HORIZON>=7,'scan horizon is at least seven')
        I=ident(2);I4=ident(4);ii=F({16:1});u=root(2)*Q(1,2)
        X=((F(),F(1)),(F(1),F()));Y=((F(),-ii),(ii,F()));Z=((F(1),F()),(F(),F(-1)))
        H=scale(((F(1),F(1)),(F(1),F(-1))),u)
        T=((F(1),F()),(F(),(1+ii)*u));S=mm(T,T);R=mm(Z,H)
        HA=kron(H,I);TA=kron(T,I);SA=kron(S,I);RB=kron(I,R)
        C=tuple(tuple(F(j==((i^1) if i>=2 else i)) for j in range(4)) for i in range(4))
        labels=[a+b for a in 'IXYZ' for b in 'IXYZ'];paulis=[kron(a,b) for a in (I,X,Y,Z) for b in (I,X,Y,Z)]
        def signed_pauli(A):
            cs=[trace(mm(P,A))*Q(1,4) for P in paulis];nz=[i for i,z in enumerate(cs) if z!=0]
            check(len(nz)==1 and cs[nz[0]] in (F(1),F(-1)),'exact signed Pauli image')
            j=nz[0];return ('-' if cs[j]==-1 else '')+labels[j]
        check(mm(adj(R),R)==I and R==scale(addmat(I,scale(Y,ii)),u),'R=ZH=(I+iY)/sqrt(2) is a unitary local Clifford')
        check(mm(mm(R,X),adj(R))==Z and mm(mm(R,Z),adj(R))==scale(X,-1),'quarter-turn Pauli action of ZH')
        delta=addmat(I,scale(R,-1))
        check(mm(adj(delta),delta)==scale(I,2-root(2)),'exact quarter-turn error square 2-sqrt(2)')
        check(trace(R)==root(2),'eigenphase pair has positive trace sqrt(2)')
        # Type-dependent Clifford multiplier. All other records stay fixed.
        multipliers={}
        for a in range(2):
            for b in range(2):
                E=kron(T if a else I,T if b else I)
                D=mm(C,RB) if (a,b)==(1,0) else C;multipliers[a,b]=D
                check(mm(adj(D),D)==I4,'each CNOT table type is a bijective Clifford left multiplication')
                for P in paulis:signed_pauli(mm(mm(D,P),adj(D)))
                V=mm(mm(E,D),adj(E));rel=mm(adj(C),V)
                if (a,b)==(1,0):
                    check(V==mm(C,RB) and rel==RB,'modified type realises CNOT after the local target quarter-turn')
                elif b==0:check(rel==I4,'calibrated CNOT phase type remains exact')
                else:
                    ds=(F(1),F(1),(1+ii)*u,(1-ii)*u)
                    expected=tuple(tuple(ds[i] if i==j else F() for j in range(4)) for i in range(4))
                    check(rel==expected,'unchanged odd-target type retains its original error spectrum')
        check(multipliers[0,0]==C and mm(C,C)==I4,'CNOT calibration I <-> C preserved')
        check(mm(C,mm(C,HA))==HA,'CNOT calibration H_A <-> C H_A preserved')
        check(mm(HA,HA)==I4,'unchanged exact Hadamard calibration')
        Hrel=mm(adj(H),mm(mm(T,H),adj(T)));hd=addmat(I,scale(Hrel,-1))
        check(mm(adj(hd),hd)==scale(I,1-u),'unchanged Hadamard error certificate')
        odd=mm(mm(C,RB),HA);even=mm(mm(C,HA),SA);U=mm(C,mm(HA,TA))
        expected_stabilisers=[['ZI','IZ'],['XX','-IX'],['-YI','-IX'],['YX','-ZZ'],['-ZX','YY'],['YY','XZ'],['-IY','-XY'],['-ZY','-IY']]
        expected_ideal=[Q(0),Q(1,2),Q(1,2),Q(1,4),Q(1,4),Q(3,8),Q(5,8),Q(1,16)]
        G=I4;ideal_matrix=I4;physical=I4;a=0;actual=[];ideal=[];joint=[];stab=[];first=None
        for n in range(SCAN_HORIZON+1):
            frame=mm(TA if a else I4,G)
            check(frame==physical,'frame carry dynamics equals direct physical elementary-gate composition')
            q=sum((sqabs(frame[j][0]) for j in (2,3)),F())
            p=sum((sqabs(ideal_matrix[j][0]) for j in (2,3)),F())
            check(set(q.c)<={0} and set(p.c)<={0},'both measured probabilities are rational')
            q=q.c.get(0,Q(0));p=p.c.get(0,Q(0));actual.append(q);ideal.append(p)
            joint.append([encode_scalar(sqabs(frame[j][0])) for j in range(4)])
            if first is None and abs(q-p)>Q(1,3):first=n
            if n<=7:
                row=[signed_pauli(mm(mm(G,paulis[j]),adj(G))) for j in (12,3)]
                check(row==expected_stabilisers[n],'hand-derived stabiliser certificate at depth '+str(n));stab.append(row)
                check(p==expected_ideal[n],'previously certified ideal probability at depth '+str(n))
            old_a=a;a=1-a
            carry=SA if old_a else I4
            G=mm(multipliers[a,0],mm(HA,mm(carry,G)))
            E=TA if a else I4;VH=mm(mm(E,HA),adj(E));VC=mm(mm(E,multipliers[a,0]),adj(E))
            physical=mm(VC,mm(VH,mm(TA,physical)))
            ideal_matrix=mm(U,ideal_matrix)
        check(actual[:8]==[Q(0)]+[Q(1,2)]*6+[Q(0)],'counterexample probability history through seven')
        check(max(abs(actual[n]-ideal[n]) for n in range(8))==Q(1,4),'entire history through seven has maximum discrepancy one quarter')
        check(abs(actual[7]-ideal[7])==Q(1,16),'seven-block discrepancy is one sixteenth')
        check(joint[7]==['1/2','1/2','0','0'],'joint seven-block output probabilities')
        check(first is None or first>7,'no strict one-third crossing through seven')
        # The record distinguishes actual verified outputs from unexecuted expectations.
        results={'counterexample_verified':True,'modified_phase_type':[1,0],
          'CNOT_Clifford_multiplier_at_modified_type':'C (I tensor ZH)','other_phase_types':'C',
          'H_error_bound':'sqrt(1-1/sqrt(2))','C_error_bound':'sqrt(2-sqrt(2))',
          'retained_labels':184320,'additional_labels_relative_to_previous_stage':0,
          'calibration_preserved':True,'inverse_C_rule':'left multiplication by D_ab^dagger; phase types unchanged',
          'actual_A1_n_0_to_7':[str(x) for x in actual[:8]],'ideal_A1_n_0_to_7':[str(x) for x in ideal[:8]],
          'maximum_discrepancy_through_seven':'1/4','discrepancy_at_seven':'1/16','joint_at_seven':joint[7],
          'signed_inner_stabilisers_n_0_to_7':stab,'scan_horizon':SCAN_HORIZON,
          'first_gap_above_one_third_within_scan':first,'actual_A1_full_scan':[str(x) for x in actual],
          'ideal_A1_full_scan':[str(x) for x in ideal],
          'all_fixed_phase_priors':'analytic consequence of inherited right-phase equivariance and compensated preparation',
          'certification_scope':'explicit full-table algebraic rule; no exhaustive optimisation over CNOT matchings and no claim about asymptotic onset'}
    except Exception:status='failed';error=traceback.format_exc()
    record={'script':source.name,'script_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'started_utc':started,'finished_utc':now(),'python':platform.python_version(),'platform':platform.platform(),'status':status,'error':error,'total_assertions':len(checks),'passed_assertions':sum(x['passed'] for x in checks),'results':results,'checks':checks}
    (folder/'cnot_variation_run_record.json').write_text(json.dumps(record,indent=2)+'\n')
    (folder/'cnot_variation_assertion_log.txt').write_text('\n'.join(('PASS: ' if c['passed'] else 'FAIL: ')+c['label'] for c in checks)+'\n'+(error or ''))
    if status=='passed':(folder/'cnot_variation_results.json').write_text(json.dumps(results,indent=2)+'\n')
    print('Mathematical verification:',status);print('Local record:',folder)
    if status=='passed':print('Maximum discrepancy through seven: 1/4; discrepancy at seven: 1/16.\nFirst strict one-third crossing within scan:',first)
    if error:print(error)
    return 0 if status=='passed' else 1
if __name__=='__main__':raise SystemExit(main())
