#!/usr/bin/env python3
"""Matched H/T controls and compensated preparation, v0.1.
Standalone Python 3.10+, standard library. Exact arithmetic in
Q(sqrt(2),sqrt(3),sqrt(5),sqrt(7),i) plus rational interval bounds.
The frozen H table was selected numerically; exact optimality is NOT claimed.
This script certifies its permutation, anchors and error <=13/50.
It saves a fresh JSON record and assertion log beside this saved script.
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

def main():
    source=Path(__file__).resolve();now=lambda:dt.datetime.now(dt.timezone.utc).isoformat()
    started=now();checks=[];summary={};status='passed';error=None
    folder=source.parent/'reproducibility'/('matched_controls_python_'+dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8]);folder.mkdir(parents=True)
    def check(ok,label):
        checks.append({'label':label,'passed':bool(ok)})
        if not ok:raise AssertionError(label)
    try:
        cat,phases=catalogue();n=len(cat);eye=ident(2);h=scale(((F(1),F(1)),(F(1),F(-1))),root(2)*Q(1,2));t=((F(1),F()),(F(),phases[1]))
        p=list(H_PERM);pt=[phase_index(i) for i in range(n)]
        check(n==464 and sorted(p)==list(range(n)),'frozen H table is a 464-label permutation')
        check(sorted(pt)==list(range(n)),'exact T table is a permutation')
        check(cat[12]==eye and cat[208]==h and p[12]==208 and p[208]==12,'exact calibrated H anchors')
        ih=[p.index(i) for i in range(n)];it=[pt.index(i) for i in range(n)]
        maps={'H':p,'h':ih,'T':pt,'t':it}
        def end(a,word):
            for g in word:a=maps[g][a]
            return a
        lows=[];threshold=(2-Q(13,50)**2)**2
        for a in range(n):
            check(mm(adj(cat[a]),cat[a])==eye,'exact catalogue unitarity')
            check(mm(t,cat[a])==cat[pt[a]],'T is exact common left multiplication')
            z=trace(mm(adj(cat[p[a]]),mm(h,cat[a])));q=sqabs(z);lo,hi=bounds(q);lows.append(lo)
            check(lo>=threshold,'rational interval certificate: matched H error at most 13/50')
            check(ih[p[a]]==a and it[pt[a]]==a,'exact inverse commands')
        words=('HTH','HTHTH','THHTHT','HTtHhT','HTHTHTHT')
        for word in words:
            inv=''.join(g.swapcase() for g in word[::-1])
            for a in range(n):check(end(end(a,word),inv)==a,'every labelled frame restored by inverse word')
        starts=[8+(k+4)%8 for k in range(8)]
        rows=[]
        for word,expected in [('HTH',[Q(7,8)]*8),('HTHTH',[Q(3,4)]*4+[Q(7,8)]*3+[Q(3,4)])]:
            probs=[];ends=[]
            for k,a in enumerate(starts):
                b=end(a,word);ends.append(b);v=mm(cat[b],adj(cat[a]));q=sqabs(v[0][0]);check(q==expected[k],'exact phase-conditioned output probability');probs.append(expected[k])
            rows.append({'word':word,'zero_based_final_indices':ends,'phase_probabilities':[str(q) for q in probs],'calibrated':str(probs[0]),'uniform_phase_average':str(sum(probs)/8)})
        ideal=eye
        for g in 'HTH':ideal=mm(h if g=='H' else t,ideal)
        check(sqabs(ideal[0][0])==(2+root(2))*Q(1,4),'exact ideal HTH probability')
        ideal=eye
        for g in 'HTHTH':ideal=mm(h if g=='H' else t,ideal)
        check(sqabs(ideal[0][0])==Q(3,4),'exact ideal HTHTH probability')
        # Direct two-sided preparation check on two Bell bases, not only a channel formula.
        imag=F({16:1});x=((F(),F(1)),(F(1),F()));y=((F(),-imag),(imag,F()));z=((F(1),F()),(F(),F(-1)))
        s=((F(1),F()),(F(),imag));rs=(eye,x,y,z);ensembles=(rs,tuple(mm(s,r) for r in rs))
        target=scale(ident(4),Q(1,4));rt=root(2)*Q(1,2)
        for wa,wb in [('HTH','HTHTH'),('HTHTH','THHTHT')]:
            averages=[scale(ident(4),0),scale(ident(4),0)]
            for ka,a in enumerate(starts):
                for kb,b in enumerate(starts):
                    ap=cat[end(a,wa)];bp=cat[end(b,wb)];va=mm(ap,adj(cat[a]));vb=mm(bp,adj(cat[b]))
                    for ei,en in enumerate(ensembles):
                        density=scale(ident(4),0)
                        for r in en:
                            psi=scale(r,rt);seed=mm(mm(adj(cat[a]),psi),transpose(adj(cat[b])))
                            actual=mm(mm(ap,seed),transpose(bp));common=mm(mm(va,psi),transpose(vb))
                            check(actual==common,'exact compensated-seed and common-channel endpoints agree')
                            density=addmat(density,scale(outer(flat(actual)),Q(1,4)))
                        check(density==target,'each calibration branch preserves both I/4 preparation recipes')
                        averages[ei]=addmat(averages[ei],scale(density,Q(1,64)))
            check(averages==[target,target],'uniform independent phase preparation is two-sided mixture consistent')
        summary={'D':8,'L':8,'N':464,'certified_H_error_bound':'13/50','minimum_trace_squared_lower_bound':str(min(lows)),
          'phase_control':'T=diag(1,exp(i*pi/4)); exact','probability_rows':rows,
          'arithmetic':'exact multiquadratic field and rational intervals; no floating point in verification',
          'scope':'Certificates for the frozen candidate; all-word feasibility follows from the general theorem, not word sampling',
          'table_sha256':hashlib.sha256(json.dumps(p,separators=(',',':')).encode()).hexdigest()}
    except Exception as exc:status='failed';error=repr(exc)
    record={'script':source.name,'script_version':'0.1','script_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'started_utc':started,'finished_utc':now(),'python_version':platform.python_version(),'platform':platform.platform(),'status':status,'error':error,'passed_assertions':sum(c['passed'] for c in checks),'total_assertions':len(checks),'summary':summary,'checks':checks}
    (folder/'matched_controls_run_record.json').write_text(json.dumps(record,indent=2)+'\n');(folder/'matched_controls_assertion_log.txt').write_text('\n'.join(('PASS: ' if c['passed'] else 'FAIL: ')+c['label'] for c in checks)+'\n')
    print(json.dumps({k:v for k,v in record.items() if k!='checks'},indent=2));print('Local record:',folder)
    if status!='passed':raise SystemExit(1)
if __name__=='__main__':main()
