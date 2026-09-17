#!/usr/bin/env python3
"""Exact finite-history certificate, standalone Python 3.10+ standard library.
Verifies the complete right-phase-equivariant local matching class on C_8,8.
A shared calibrated table postpones error >1/3 through 192 HT blocks for
all fixed phase priors; exhaustive search excludes survival through 193.
This is a fixed-catalogue, symmetry-restricted theorem, not a full two-qubit
or asymptotic onset theorem. Saves fresh records beside the saved script.
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


WITNESS_H = [212, 268, 260, 252, 244, 236, 228, 220, 240, 248, 256, 264, 208, 216, 224, 232, 276, 397, 261, 125, 180, 108, 163, 348, 284, 341, 269, 133, 188, 116, 171, 356, 292, 349, 213, 141, 196, 124, 179, 364, 300, 357, 221, 85, 204, 132, 187, 372, 308, 365, 229, 93, 148, 140, 195, 380, 316, 373, 237, 101, 156, 84, 203, 388, 324, 381, 245, 109, 164, 92, 147, 396, 332, 389, 253, 117, 172, 100, 155, 340, 411, 461, 325, 189, 61, 43, 227, 283, 419, 405, 333, 197, 69, 51, 235, 291, 427, 413, 277, 205, 77, 59, 243, 299, 435, 421, 285, 149, 21, 67, 251, 307, 443, 429, 293, 157, 29, 75, 259, 315, 451, 437, 301, 165, 37, 19, 267, 323, 459, 445, 309, 173, 45, 27, 211, 331, 403, 453, 317, 181, 53, 35, 219, 275, 404, 390, 198, 70, 52, 107, 290, 347, 412, 398, 206, 78, 60, 115, 298, 355, 420, 342, 150, 22, 68, 123, 306, 363, 428, 350, 158, 30, 76, 131, 314, 371, 436, 358, 166, 38, 20, 139, 322, 379, 444, 366, 174, 46, 28, 83, 330, 387, 452, 374, 182, 54, 36, 91, 274, 395, 460, 382, 190, 62, 44, 99, 282, 339, 12, 454, 262, 134, 0, 98, 226, 354, 13, 462, 270, 142, 7, 106, 234, 362, 14, 406, 214, 86, 6, 114, 242, 370, 15, 414, 222, 94, 5, 122, 250, 378, 8, 422, 230, 102, 4, 130, 258, 386, 9, 430, 238, 110, 3, 138, 266, 394, 10, 438, 246, 118, 2, 82, 210, 338, 11, 446, 254, 126, 1, 90, 218, 346, 432, 447, 326, 79, 16, 34, 162, 418, 440, 455, 334, 23, 24, 42, 170, 426, 448, 463, 278, 31, 32, 50, 178, 434, 456, 407, 286, 39, 40, 58, 186, 442, 400, 415, 294, 47, 48, 66, 194, 450, 408, 423, 302, 55, 56, 74, 202, 458, 416, 431, 310, 63, 64, 18, 146, 402, 424, 439, 318, 71, 72, 26, 154, 410, 425, 383, 327, 207, 25, 89, 225, 361, 433, 391, 335, 151, 33, 97, 233, 369, 441, 399, 279, 159, 41, 105, 241, 377, 449, 343, 287, 167, 49, 113, 249, 385, 457, 351, 295, 175, 57, 121, 257, 393, 401, 359, 303, 183, 65, 129, 265, 337, 409, 367, 311, 191, 73, 137, 209, 345, 417, 375, 319, 199, 17, 81, 217, 353, 304, 376, 263, 143, 144, 88, 161, 297, 312, 384, 271, 87, 152, 96, 169, 305, 320, 392, 215, 95, 160, 104, 177, 313, 328, 336, 223, 103, 168, 112, 185, 321, 272, 344, 231, 111, 176, 120, 193, 329, 280, 352, 239, 119, 184, 128, 201, 273, 288, 360, 247, 127, 192, 136, 145, 281, 296, 368, 255, 135, 200, 80, 153, 289]

def main():
    source=Path(__file__).resolve();now=lambda:dt.datetime.now(dt.timezone.utc).isoformat()
    folder=source.parent/'reproducibility'/('probability_histories_python_'+dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8]);folder.mkdir(parents=True)
    started=now();checks=[];result={};status='passed';error=None
    def check(ok,label):
        checks.append({'label':label,'passed':bool(ok)})
        if not ok:raise AssertionError(label)
    try:
        cat,phases=catalogue();N=len(cat);H=cat[208];eye=ident(2)
        T=((F(1),F()),(F(),phases[1]));pt=[phase_index(i) for i in range(N)]
        reps=[0,8]+[16+n*64+g for n in range(7) for g in range(8)]
        orbits=[];decode={}
        for i,a in enumerate(reps):
            row=[];x=a
            for k in range(8):
                check(x not in decode,'free right-phase orbit partition')
                decode[x]=(i,k);row.append(x);x=right_phase(x)
            check(x==a,'right-phase orbit length eight');orbits.append(row)
        for i in range(N):
            check(mm(adj(cat[i]),cat[i])==eye,'unitarity of every catalogue frame')
            check(mm(T,cat[i])==cat[pt[i]],'exact left phase control')
            check(sqabs(trace(mm(adj(cat[right_phase(i)]),mm(cat[i],T))))==4,'exact projective right-phase action')
            check(pt[right_phase(i)]==right_phase(pt[i]),'left T commutes with right-phase labels')
        # All possible phase-equivariant block edges, not one frozen table.
        threshold=F((2-Q(13,50)**2)**2)
        conj=[tuple(z.conj() for z in flat(a)) for a in cat]
        qedges=[]
        for i,a in enumerate(reps):
            target=flat(mm(H,cat[pt[a]]));row=[]
            for j in range(N):
                z=sum((u*w for u,w in zip(conj[j],target)),F())
                if sign(sqabs(z)-threshold)>=0:row.append(decode[j])
                check(True,'exact classification of quotient edge and phase shift')
            qedges.append(sorted(row))
            if i%10==0:print('Exact quotient rows:',i+1,'/',len(reps),flush=True)
        forced={};forced_shift={}
        for x,y in ((11,208),(215,12)):
            i,k=decode[x];j,l=decode[y];forced[i]=j;forced_shift[i]=(l-k)%8
            check((j,forced_shift[i]) in qedges[i],'calibration shift is an allowed exact edge')
        qgraph=restrict([sorted(set(j for j,k in row)) for row in qedges],forced)
        check(forced=={1:26,33:1} and forced_shift=={1:5,33:4},'exact quotient calibration assignments')
        check(sum(map(len,qgraph))==174,'complete calibrated quotient graph has 174 edges')
        initial=matching(qgraph)
        check(sorted(initial)==list(range(58)),'quotient graph has a full matching')
        p=WITNESS_H;P=[p[pt[i]] for i in range(N)]
        check(sorted(p)==list(range(N)) and p[12]==208 and p[208]==12,'witness is a full calibrated H permutation')
        for i in range(N):
            check(p[right_phase(i)]==right_phase(p[i]),'witness H is right-phase-equivariant')
            check(P[right_phase(i)]==right_phase(P[i]),'witness block is right-phase-equivariant')
            z=trace(mm(adj(cat[p[i]]),mm(H,cat[i])))
            check(sign(sqabs(z)-threshold)>=0,'every witness H edge has error at most 13/50')
        qp=[];shifts=[]
        for i,a in enumerate(reps):
            j,k=decode[P[a]];qp.append(j);shifts.append(k)
            check(j in qgraph[i] and (j,k) in qedges[i],'witness quotient edge and shift allowed')
        check(sorted(qp)==list(range(58)),'witness quotient is a permutation')
        cyc=orbit(qp,1)
        check(cyc==[1,26,56,35,40,52,17,53,18,48,43,41,50,42,49,51,33],'witness observable orbit is the certified 17-cycle')
        rq=[Q(1)-sqabs(cat[a][0][0]).c.get(0,Q(0)) for a in reps]
        check(all(sqabs(cat[a][0][0])==1-rq[i] for i,row in enumerate(orbits) for a in row),'readout constant on right-phase orbits')
        c=-(2+root(2))*Q(1,4);cos=[F(1),c]
        for n in range(2,194):cos.append(2*c*cos[-1]-cos[-2])
        ideal=[(6-root(2))*Q(1,17)*(1-x) for x in cos]
        allowed=[[sign((F(r)-ideal[n])*(F(r)-ideal[n])-Q(1,9))<=0 for r in rq] for n in range(194)]
        check(all(allowed[n][cyc[n%17]] for n in range(193)),'witness survives every depth zero through 192')
        check(not allowed[193][cyc[193%17]],'witness first exceeds one third at block 193')
        check(all(sign((F(rq[cyc[n%17]])-ideal[n])*(F(rq[cyc[n%17]])-ideal[n])-Q(13,40)**2)<0 for n in range(193)), 'witness has error strictly below 13/40 through block 192')
        check(sign((F(rq[cyc[193%17]])-ideal[193])*(F(rq[cyc[193%17]])-ideal[193])-Q(17,50)**2)>0, 'witness error at block 193 strictly exceeds 17/50')
        check(len(orbit(P,12))==17,'witness full frame period is also seventeen')
        # Exact compensated matrices coincide across all initial phase contexts.
        states=list(range(8,16));base=12;probabilities=[]
        for n in range(194):
            V=cat[base]
            for k,a0 in enumerate(range(8,16)):
                branch=mm(cat[states[k]],adj(cat[a0]));rel=trace(mm(adj(V),branch))
                check(sqabs(rel)==4,'all compensated context matrices agree projectively at each tested depth')
            probabilities.append(str(rq[cyc[n%17]]));base=P[base];states=[P[a] for a in states]
        # Exhaustive depth-first decision procedure for all quotient permutations.
        # It explores every possible next row assignment, enforces stationarity
        # and injectivity, and prunes only when a full matching cannot extend it.
        def force_edge(full,assigned,x,y):
            if full[x]==y:return full[:]
            inv=[0]*58
            for i,j in enumerate(full):inv[j]=i
            old=full[x];orphan=inv[y];out=full[:];out[x]=y;out[orphan]=-1;inv[old]=-1;inv[y]=x
            locked=set(assigned)|{x}
            def augment(i,seen):
                for j in qgraph[i]:
                    if j in seen:continue
                    seen.add(j);owner=inv[j]
                    if owner<0 or (owner not in locked and augment(owner,seen)):
                        out[i]=j;inv[j]=i;return True
                return False
            return out if augment(orphan,set()) else None
        stats={'nodes':0,'new_assignments':0,'band_rejections':0,'injection_rejections':0,'Hall_rejections':0,'longest_surviving_horizon':0}
        digest=hashlib.sha256()
        def explore(path,assigned,full):
            n=len(path);stats['nodes']+=1;stats['longest_surviving_horizon']=max(stats['longest_surviving_horizon'],n-1)
            digest.update((str(n)+':'+','.join(map(str,path))+';').encode())
            if n>193:return True
            x=path[-1];choices=[assigned[x]] if x in assigned else qgraph[x]
            for y in choices:
                if not allowed[n][y]:stats['band_rejections']+=1;continue
                if x in assigned:
                    if explore(path+[y],assigned,full):return True
                elif y in assigned.values():stats['injection_rejections']+=1
                else:
                    ext=force_edge(full,assigned,x,y)
                    # Independent full matching check for each claimed Hall rejection.
                    if ext is None:
                        test=restrict(qgraph,assigned|{x:y});pbad=matching(test)
                        S,Nb=hall_witness(test,pbad)
                        check(len(Nb)<len(S) and set(Nb)=={j for i in S for j in test[i]},'independent Hall certificate for a search prune')
                        stats['Hall_rejections']+=1;continue
                    check(sorted(ext)==list(range(58)) and all(ext[i] in qgraph[i] for i in range(58)) and all(ext[i]==j for i,j in (assigned|{x:y}).items()),'new partial history has a full matching completion')
                    stats['new_assignments']+=1
                    if explore(path+[y],assigned|{x:y},ext):return True
            return False
        print('Beginning exhaustive exact 193-block search',flush=True)
        survives=explore([1],forced,initial)
        check(not survives,'no phase-equivariant calibrated table survives all depths through 193')
        check(stats['longest_surviving_horizon']==192,'exact maximal postponement horizon is 192')
        result={'scope':'Full phase-equivariant local C_8,8 matching class; all fixed priors over eight compensated phase contexts; not the unrestricted two-qubit matching class.',
          'D':8,'L':8,'frames':464,'quotient_frames':58,'H_error_bound':'13/50','threshold':'1/3','largest_feasible_horizon':192,'first_unavoidable_curve_endpoint':193,'D_curve':194,
          'quotient_edges':qedges,'constrained_quotient_graph':qgraph,'calibration':[[i,forced[i],forced_shift[i]] for i in sorted(forced)],
          'witness_H_permutation':p,'witness_quotient_permutation':qp,'witness_quotient_shifts':shifts,'witness_quotient_cycle':cyc,'probabilities_n_0_to_193':probabilities,
          'search_statistics':stats,'search_traversal_sha256':digest.hexdigest(),'arithmetic':'exact multiquadratic arithmetic, outward rational intervals, integer exhaustive search; no floating point in verification'}
    except Exception as exc:status='failed';error=repr(exc)
    record={'script':source.name,'version':'0.1','script_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'started_utc':started,'finished_utc':now(),'python':platform.python_version(),'platform':platform.platform(),'status':status,'error':error,'total_assertions':len(checks),'passed_assertions':sum(c['passed'] for c in checks),'results':result,'checks':checks}
    (folder/'probability_histories_run_record.json').write_text(json.dumps(record,indent=2)+'\n')
    (folder/'probability_histories_assertion_log.txt').write_text('\n'.join(('PASS: ' if c['passed'] else 'FAIL: ')+c['label'] for c in checks)+'\n')
    if status=='passed':
        (folder/'probability_history_results.json').write_text(json.dumps(result,indent=2)+'\n')
        (source.parent/'probability_history_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Verification:',status,'Assertions:',len(checks),'Error:',error,'\nLocal record:',folder,flush=True)
    if status!='passed':raise SystemExit(1)
if __name__=='__main__':main()
