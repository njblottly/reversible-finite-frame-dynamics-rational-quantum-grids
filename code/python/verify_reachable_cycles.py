#!/usr/bin/env python3
"""Reachable-cycle certificates. Python 3.10+, standard library only.
Exact multiquadratic arithmetic and outward rational intervals; no floats.
Standalone: run this saved file. Fresh records are saved beside it.
"""
from fractions import Fraction as Q
from pathlib import Path
from math import isqrt
import datetime as dt
import hashlib,json,platform,uuid

# Replaced by the frozen zero-based table when this script is authored.
H_PERM = (212, 268, 260, 252, 244, 236, 228, 220, 240, 248, 256, 264, 208, 216, 224, 232, 276, 396, 197, 124, 180, 43, 163, 283, 284, 277, 333, 132, 188, 51, 299, 291, 292, 348, 149, 140, 196, 59, 307, 427, 300, 356, 285, 84, 204, 67, 187, 435, 308, 364, 165, 92, 148, 75, 323, 315, 316, 372, 301, 100, 156, 19, 203, 388, 324, 380, 309, 108, 164, 27, 275, 331, 332, 325, 317, 116, 172, 35, 155, 340, 461, 397, 261, 189, 61, 171, 227, 347, 405, 341, 269, 133, 69, 179, 235, 355, 413, 349, 213, 205, 77, 123, 243, 363, 421, 293, 221, 85, 21, 195, 251, 371, 429, 365, 229, 157, 29, 139, 259, 379, 437, 373, 237, 101, 37, 147, 267, 387, 445, 381, 245, 173, 45, 91, 211, 395, 453, 389, 253, 181, 53, 99, 219, 339, 404, 390, 198, 125, 52, 107, 162, 411, 412, 398, 206, 78, 60, 115, 170, 419, 420, 342, 150, 141, 68, 50, 178, 370, 428, 357, 158, 30, 76, 131, 186, 378, 436, 358, 166, 93, 20, 66, 194, 443, 444, 366, 174, 46, 28, 83, 202, 451, 452, 374, 182, 109, 36, 18, 146, 459, 460, 382, 190, 117, 44, 26, 154, 403, 12, 454, 262, 134, 0, 98, 226, 418, 13, 462, 270, 142, 7, 106, 234, 426, 14, 406, 214, 86, 6, 114, 242, 434, 15, 414, 222, 94, 5, 122, 250, 442, 8, 422, 230, 102, 4, 130, 258, 450, 9, 430, 238, 110, 3, 138, 266, 458, 10, 438, 246, 118, 2, 82, 210, 402, 11, 446, 254, 126, 1, 90, 218, 410, 432, 383, 326, 70, 16, 34, 290, 354, 440, 391, 334, 23, 24, 42, 298, 362, 448, 399, 278, 22, 32, 41, 306, 377, 456, 350, 286, 39, 40, 58, 314, 385, 400, 351, 294, 38, 48, 57, 322, 386, 408, 359, 302, 55, 56, 74, 330, 394, 416, 367, 310, 54, 64, 73, 274, 338, 424, 375, 318, 62, 72, 17, 282, 346, 425, 447, 327, 143, 79, 89, 161, 361, 433, 455, 207, 87, 33, 97, 297, 369, 463, 335, 151, 95, 31, 105, 305, 313, 407, 279, 159, 103, 49, 113, 185, 321, 415, 287, 295, 111, 47, 121, 193, 393, 401, 423, 175, 119, 65, 129, 329, 337, 409, 303, 183, 127, 63, 137, 273, 281, 417, 311, 319, 135, 71, 81, 153, 289, 304, 376, 263, 136, 144, 25, 225, 360, 312, 384, 271, 80, 152, 96, 233, 368, 320, 392, 215, 88, 160, 169, 241, 441, 328, 343, 223, 167, 168, 177, 249, 449, 272, 344, 231, 104, 176, 120, 257, 457, 280, 352, 239, 112, 184, 128, 265, 336, 288, 431, 247, 191, 192, 201, 209, 345, 296, 439, 255, 199, 200, 145, 217, 353)

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

def bounds(x,bits=80):
    assert all(not k&16 for k in x.c),'real field element required'
    lo=hi=Q(0);m=2**bits
    for mask,c in x.c.items():
        n=FACTORS[mask];k=isqrt(n*m*m);a=Q(k,m);b=a if k*k==n*m*m else Q(k+1,m)
        lo+=c*(a if c>=0 else b);hi+=c*(b if c>=0 else a)
    return lo,hi

def catalogue():
    ii=F({16:1});u=root(2)*Q(1,2)
    phases=(F(1),u*(1+ii),ii,u*(-1+ii),F(-1),u*(-1-ii),-ii,u*(1-ii))
    def v(n,b,g):
        a=root(2*n)*Q(1,4);c=root(2*(8-n))*Q(1,4)
        return ((a,c*phases[b%8]),(c*phases[g%8],-a*phases[(b+g)%8]))
    cat=[v(0,0,k) for k in range(8)]+[v(8,0,k) for k in range(8)]
    cat += [v(n,b,g) for n in range(1,8) for b in range(8) for g in range(8)]
    return cat,phases

def phase_index(i):
    return (i//8)*8+(i%8+1)%8  # gamma shifts by one for all catalogue blocks

from collections import deque

def sign(x):
    if x==0:return 0
    for bits in (80,160,320,640,1280,2560):
        lo,hi=bounds(x,bits)
        if lo>0:return 1
        if hi<0:return -1
    raise ArithmeticError('Sign unresolved; increase certified precision')

def right_phase(i):
    if i<8:return (i-1)%8
    if i<16:return 8+(i-8+1)%8
    a=i-16;n,b,g=a//64,(a%64)//8,a%8
    return 16+n*64+((b+1)%8)*8+g

def scc(graph):
    n=len(graph);rev=[[] for _ in graph]
    for i,row in enumerate(graph):
        for j in row:rev[j].append(i)
    seen=set();order=[]
    def visit(x):
        seen.add(x)
        for y in graph[x]:
            if y not in seen:visit(y)
        order.append(x)
    for x in range(n):
        if x not in seen:visit(x)
    labels=[-1]*n;label=0
    for x in reversed(order):
        if labels[x]>=0:continue
        stack=[x];labels[x]=label
        while stack:
            u=stack.pop()
            for v in rev[u]:
                if labels[v]<0:labels[v]=label;stack.append(v)
        label+=1
    return labels

def restrict(graph,forced):
    if len(set(forced.values()))!=len(forced):return None
    if any(j not in graph[i] for i,j in forced.items()):return None
    cols=set(forced.values())
    return [[forced[i]] if i in forced else [j for j in row if j not in cols]
            for i,row in enumerate(graph)]

def matching(graph):
    n=len(graph);owner=[-1]*n
    def aug(i,seen):
        for j in graph[i]:
            if j in seen:continue
            seen.add(j)
            if owner[j]<0 or aug(owner[j],seen):owner[j]=i;return True
        return False
    for i in range(n):aug(i,set())
    p=[-1]*n
    for j,i in enumerate(owner):
        if i>=0:p[i]=j
    return p

def hall_witness(graph,p):
    owner={j:i for i,j in enumerate(p) if j>=0}
    left=set(i for i,j in enumerate(p) if j<0);todo=list(left);right=set()
    while todo:
        i=todo.pop()
        for j in graph[i]:
            if j in right:continue
            right.add(j)
            if j in owner and owner[j] not in left:
                left.add(owner[j]);todo.append(owner[j])
    return sorted(left),sorted(right)

def orbit(p,start):
    out=[];x=start
    while x not in out:out.append(x);x=p[x]
    assert x==start
    return out

def main():
    source=Path(__file__).resolve();started=dt.datetime.now(dt.timezone.utc).isoformat()
    folder=source.parent/'reproducibility'/('reachable_cycles_python_'+dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8]);folder.mkdir(parents=True)
    checks=[];result={};status='passed';error=None
    def check(ok,label):
        checks.append({'label':label,'passed':bool(ok)})
        if not ok:raise AssertionError(label)
    try:
        cat,phases=catalogue();N=len(cat);eye=ident(2)
        H=cat[208];T=((F(1),F()),(F(),phases[1]));threshold=F((2-Q(13,50)**2)**2)
        check(cat[12]==eye,'identity label and calibrated Hadamard')
        hc=[flat(mm(H,a)) for a in cat];conj=[tuple(z.conj() for z in flat(a)) for a in cat]
        for i in range(N):
            check(mm(adj(cat[i]),cat[i])==eye,'exact unitary catalogue')
            check(mm(T,cat[i])==cat[phase_index(i)],'exact left phase action')
            z=trace(mm(adj(cat[right_phase(i)]),mm(cat[i],T)))
            check(sqabs(z)==4,'right phase symmetry up to global phase')
        reps=[0,8]+[16+n*64+g for n in range(7) for g in range(8)]
        graph=[[] for _ in cat];assigned=set()
        for ri,rep in enumerate(reps):
            row=[]
            for j in range(N):
                z=sum((a*b for a,b in zip(conj[j],hc[rep])),F())
                good=sign(sqabs(z)-threshold)>=0
                if good:row.append(j)
                check(True,'exact representative edge classification')
            i=rep
            for k in range(8):
                check(i not in assigned,'right symmetry orbits partition catalogue')
                assigned.add(i);graph[i]=sorted(row)
                i=right_phase(i);row=[right_phase(j) for j in row]
            if ri%10==0:print('Exact graph representatives:',ri+1,'/',len(reps),flush=True)
        graph=restrict(graph,{12:208,208:12})
        check(graph is not None,'calibration anchors compatible with geometry')
        p=list(H_PERM)
        check(sorted(p)==list(range(N)) and all(p[i] in graph[i] for i in range(N)),'existing full matching certifies feasibility')
        inverse=[p.index(j) for j in range(N)]
        labels=scc([[inverse[j] for j in row] for row in graph])
        supported=[[j for j in row if labels[i]==labels[inverse[j]]] for i,row in enumerate(graph)]
        check(sum(map(len,graph))==1420,'calibrated geometric graph: 1420 edges')
        check(sum(map(len,supported))==1264,'supported graph: 1264 edges')
        unsupported=next((i,j) for i,row in enumerate(graph) for j in row if j not in supported[i])
        bad=restrict(graph,{unsupported[0]:unsupported[1]});badp=matching(bad)
        left,right=hall_witness(bad,badp)
        check(len(right)<len(left) and set(right)=={j for i in left for j in bad[i]},'explicit Hall obstruction for unsupported edge')
        bg=[supported[phase_index(i)] for i in range(N)];bc=scc(bg)
        check(len(set(bc))==1,'supported block graph is one 464-vertex SCC')
        dist={12:0};todo=deque([12]);short=None
        while todo:
            i=todo.popleft()
            for j in bg[i]:
                if j==12:short=dist[i]+1;todo.clear();break
                if j not in dist:dist[j]=dist[i]+1;todo.append(j)
        check(short==5,'no supported block cycle through calibrated label shorter than five')
        path=[12,216,398,353,215]
        forced={12:208,208:12}
        for i,j in zip(path,path[1:]+[12]):forced[phase_index(i)]=j
        restricted=restrict(graph,forced);newh=matching(restricted)
        check(sorted(newh)==list(range(N)) and all(newh[i] in graph[i] for i in range(N)),'five-cycle extends to full calibrated H permutation')
        oldb=[p[phase_index(i)] for i in range(N)];newb=[newh[phase_index(i)] for i in range(N)]
        check(orbit(newb,12)==path and len(orbit(oldb,12))==312,'exact cycle lengths 5 and 312')
        # cos(2 n theta_s), with sin^2(theta_s)=(6+s sqrt(2))/8.
        def cosseq(s,count):
            c=-(2+s*root(2))*Q(1,4);seq=[F(1),c]
            for n in range(2,count+1):seq.append(2*c*seq[-1]-seq[-2])
            return seq
        cp=cosseq(1,512);cm=cosseq(-1,512)
        scans=[]
        for den,oldcount,oldfirst,newcount,newfirst in [(16,245,7,241,15),(64,216,11,192,22),(256,95,11,20,128),(1024,22,95,1,172)]:
            old=[];strong=[]
            for q in range(1,257):
                e2=Q(q*q,den*den);sp=(1-cp[q])*Q(1,2);sm=(1-cm[q])*Q(1,2)
                single=sign(F(e2)-sm)>=0
                # cos(q theta+)cos(q theta-) term only needed at even q.
                mu=2*(2+cp[q]+cm[q])
                if q%2==0:mu+=8*((-1)**(q//2))*cp[q//2]*cm[q//2]
                check(sign(mu)>=0,'nonnegative exact trace squared')
                full=single and sign(F(e2)-sp)>=0 and (e2>=2 or sign(mu-(4-2*e2)**2)>=0)
                if single:old.append(q)
                if full:strong.append(q)
            check((len(old),old[0],len(strong),strong[0])==(oldcount,oldfirst,newcount,newfirst),'exact spectral-screen count and first survivor')
            scans.append({'error':str(Q(1,den)),'single_count':len(old),'single_first':old[0],'strong_count':len(strong),'strong_first':strong[0],'strong_survivors':strong})
        rows=[]
        for name,b in [('frozen',oldb),('five_cycle',newb)]:
            idx=12;first={Q(1,20):None,Q(1,3):None};probs=[]
            for n in range(65):
                actual=1-sqabs(cat[idx][0][0]);ideal=(6-root(2))*Q(1,17)*(1-cp[n]);diff=actual-ideal
                # 2/(6+sqrt2)=(6-sqrt2)/17.
                for delta in first:
                    if first[delta] is None and sign(diff*diff-delta*delta)>0:first[delta]=n
                probs.append(str(actual.c.get(0,0)));idx=b[idx]
            expected=(4,6) if name=='frozen' else (2,19)
            check(tuple(first.values())==expected,'exact first strict discrepancy thresholds: '+name)
            rows.append({'table':name,'calibrated_cycle_length':len(orbit(b,12)),'first_above_1_over_20':first[Q(1,20)],'first_above_1_over_3':first[Q(1,3)],'probabilities_n_0_to_64':probs})
        result={'arithmetic':'exact multiquadratic field with outward rational intervals; zero-based labels','N':N,'D':8,'L':8,'H_error_bound':'13/50','geometric_edges_after_calibration':1420,'supported_edges':1264,'block_scc_sizes':[464],
          'unsupported_edge':list(unsupported),'forced_edge_Hall_witness':{'left':left,'neighbors':right},'shortest_calibrated_block_cycle':path,'five_cycle_H_permutation':newh,'frozen_H_permutation':p,'spectral_screens':scans,'onset_rows':rows,
          'scope':'Local (T,H) catalogue benchmark; no claim of an asymptotic or full two-qubit discrepancy onset. Spectral screens are necessary, not sufficient.'}
    except Exception as exc:status='failed';error=repr(exc)
    record={'script':source.name,'version':'0.1','script_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'started_utc':started,'finished_utc':dt.datetime.now(dt.timezone.utc).isoformat(),'python':platform.python_version(),'platform':platform.platform(),'status':status,'error':error,'total_assertions':len(checks),'passed_assertions':sum(c['passed'] for c in checks),'results':result,'checks':checks}
    (folder/'reachable_cycles_run_record.json').write_text(json.dumps(record,indent=2)+'\n')
    (folder/'reachable_cycles_assertion_log.txt').write_text('\n'.join(('PASS: ' if c['passed'] else 'FAIL: ')+c['label'] for c in checks)+'\n')
    if status=='passed':
        (folder/'reachable_cycle_results.json').write_text(json.dumps(result,indent=2)+'\n')
        (source.parent/'reachable_cycle_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Verification:',status,'Assertions:',len(checks),'Error:',error,'\nLocal record:',folder,flush=True)
    if status!='passed':raise SystemExit(1)
if __name__=='__main__':main()
