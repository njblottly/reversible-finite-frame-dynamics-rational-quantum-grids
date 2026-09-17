#!/usr/bin/env python3
"""Exact CNOT/local-frame benchmark v0.1; Python 3.10+, standard library only.
Run: python3 verify_entangling_controls.py
This exhausts the two-qubit Clifford controller, NOT a refined H/T/CNOT table.
The general refining construction is proved in raqm_extension.tex.
"""
from collections import deque
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
import hashlib, json, platform, sys, traceback, uuid

VERSION = '0.1'
ROOT = Path(__file__).resolve().parent
CHECKS = []
SUMMARY = {}
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

def main():
    # Check preservation of all 16x16 Pauli products, including signs.
    for name,p in GENS.items():
        check(sorted(map(abs,p))==list(ID) and p[0]==1, name+' signed Pauli permutation')
        valid=True
        for a in range(16):
            for b in range(16):
                z,c=mul2(a,b); x,y=p[a],p[b];zz,cc=mul2(abs(x)-1,abs(y)-1)
                zz=(zz+(2 if x<0 else 0)+(2 if y<0 else 0))%4
                rhs=p[c];z=(z+(2 if rhs<0 else 0))%4
                valid &= (cc==abs(rhs)-1 and zz==z)
        check(valid,name+' preserves every Pauli product')
    check(compose(C,C)==ID,'CNOT is an involution')
    check(C[4]==6 and C[3]==16,'CNOT maps XI to XX and IZ to ZZ')
    group,index=enumerate_group(GENS.values())
    loc,locindex=enumerate_group((HA,HB,SA,SB))
    check(len(group)==11520,'complete generated Clifford group has 11520 elements')
    check(len(loc)==576,'local subgroup has 576 elements')
    reps=[]; coordinates={}
    for g in group:
        if g in coordinates:continue
        j=len(reps);reps.append(g)
        for k,a in enumerate(loc):
            frame=compose(a,g)
            check(frame not in coordinates,'left-coset factorisation is unique')
            coordinates[frame]=(k,j)
    check(len(reps)==20 and len(coordinates)==11520,'20 left cosets cover the controller')
    check(reps[0]==ID,'preparation uses identity chart')
    tables={}
    for name,p in GENS.items():
        table=[index[compose(p,g)] for g in group];tables[name]=table
        check(sorted(table)==list(range(len(group))),name+' is a permutation on every frame')
        check(all(compose(inv(p),group[table[i]])==g for i,g in enumerate(group)),name+' exact inverse on every frame')
        if name!='C':
            check(all(coordinates[group[table[i]]][1]==coordinates[g][1] for i,g in enumerate(group)),name+' preserves nonlocal chart on every frame')
            check(all(compose(group[table[i]],inv(g))==p for i,g in enumerate(group)),name+' physical update is the same strictly local gate')
    check(all(tables['C'][tables['C'][i]]==i for i in range(11520)),'CNOT echo restores every frame')
    check(all(compose(HA,compose(HB,g))==compose(HB,compose(HA,g)) for g in group),'opposite local Hadamards commute on every frame')
    check(C not in locindex and compose(C,compose(HA,C)) not in locindex,'CNOT and conjugated local Hadamard are nonlocal')
    # Density coefficients r_j=Tr(rho P_j); density matrix is sum_j r_j P_j/4.
    z00=state({3:1,12:1,15:1}); phi=state({5:1,10:-1,15:1})
    plusplus=state({1:1,4:1,5:1}); plusi0=state({3:1,8:1,11:1})
    check(act(compose(C,HA),z00)==phi,'calibrated H_A then CNOT creates Bell state exactly')
    check(all(phi[j]==0 for j in (1,2,3,4,8,12)),'Bell marginals are exactly maximally mixed')
    for a in range(2):
        for b in range(2):
            inp=state({12:(-1)**a,3:(-1)**b,15:(-1)**(a+b)})
            out=state({12:(-1)**a,3:(-1)**(a^b),15:(-1)**(a+(a^b))})
            check(act(C,inp)==out,'CNOT computational truth table')
    # Two Bell decompositions of I/4, and several coherent pure inputs.
    paulis_a=[ID,local((1,2,-3,-4),0),local((1,-2,3,-4),0),local((1,-2,-3,4),0)]
    bell=[act(p,phi) for p in paulis_a]
    rotated=[act(SA,r) for r in bell]
    mixed=state({})
    check(average(bell)==mixed and average(rotated)==mixed,'two preparation recipes have equal input density')
    words=[['HA','C'],['HA','C','SB','HB','C','SA'],['HB','SA','C','HA','SB','C','HB']]
    phases=[compose(power(SA,k),power(SB,l)) for k in range(4) for l in range(4)]
    probability_rows=[]
    for word in words:
        u=word_product(word)
        check(all(compose(inv(u),compose(u,g))==g for g in group),'entangling word inverse restores every frame: '+' '.join(word))
        averages=[]
        for f0 in phases:
            frame=compose(u,f0);w=compose(frame,inv(f0))
            check(w==u,'compensated common unitary independent of preparation phase')
            for r in [z00,phi,plusplus,plusi0]+bell+rotated:
                seed=act(inv(f0),r)
                direct=act(frame,seed)
                check(direct==act(w,r),'direct compensated seed equals common-channel output')
                check(act(inv(frame),direct)==seed,'seed remains recoverable after entangling word')
            for recipe in (bell,rotated):
                averages.append(average([act(frame,act(inv(f0),r)) for r in recipe]))
        check(all(r==mixed for r in averages),'both I/4 recipes remain equal in every phase branch')
        rr=act(u,z00)
        probs=[Fraction(1+(-1)**a*rr[12]+(-1)**b*rr[3]+(-1)**(a+b)*rr[15],4) for a in range(2) for b in range(2)]
        check(sum(probs)==1 and min(probs)>=0,'normalised exact computational readout')
        probability_rows.append({'chronological_word':word,'probabilities':[str(p) for p in probs]})
    SUMMARY.update({'group_size':len(group),'local_subgroup_size':len(loc),'left_cosets':len(reps),
        'D':2,'L':4,'single_qubit_frames':24,'initial_phase_pairs':16,
        'complete_generator_transition_checks':5*len(group),'probability_rows':probability_rows,
        'scope':'Exact Clifford benchmark only. No non-Clifford CNOT matching table or numerical quotient-radius constant certified.',
        'general_result':'The refining H/T/CNOT construction and all-word preparation law are analytic theorems in the manuscript.'})

if __name__=='__main__':
    started=datetime.now(timezone.utc)
    out=ROOT/'reproducibility'/('entangling_python_'+started.strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8])
    out.mkdir(parents=True)
    status='failed';error=None
    try:
        main();status='passed'
    except Exception:
        error=traceback.format_exc();print(error)
    record={'script':Path(__file__).name,'script_version':VERSION,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'started_utc':started.isoformat(),'finished_utc':datetime.now(timezone.utc).isoformat(),
        'python_version':sys.version,'platform':platform.platform(),'arithmetic':'exact integers and fractions; no floating point',
        'status':status,'passed_assertions':sum(c['passed'] for c in CHECKS),'total_assertions':len(CHECKS),
        'summary':SUMMARY,'checks':CHECKS,'error':error}
    (out/'entangling_controls_run_record.json').write_text(json.dumps(record,indent=2)+'\n')
    (out/'entangling_controls_assertion_log.txt').write_text('\n'.join(('PASS: ' if c['passed'] else 'FAIL: ')+c['label'] for c in CHECKS)+'\n'+json.dumps(SUMMARY,indent=2)+'\n'+(error or ''))
    print(json.dumps(SUMMARY,indent=2));print('Mathematical verification:',status)
    print('Assertions:',record['passed_assertions'],'/',len(CHECKS));print('Local record:',out)
    sys.exit(0 if status=='passed' else 1)
