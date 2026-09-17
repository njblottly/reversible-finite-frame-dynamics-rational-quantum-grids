#!/usr/bin/env python3
"""Checks for Nearby Admissible Gates and Precision Allocation.

Python 3.10+, standard library only. Exact fraction checks are explicitly
separated from numerical matrix checks (tolerance 1e-10). This program is
not an interval certificate and does not simulate physical values of L.
"""
from __future__ import annotations
import argparse
import cmath
import json
import math
import random
from fractions import Fraction as F
from itertools import product
from pathlib import Path

TOL = 1e-10
SEED = 20260907

def inner(x, y):
    return sum(a.conjugate()*b for a,b in zip(x,y))

def norm(x):
    return math.sqrt(sum(abs(a)**2 for a in x))

def normalise(x):
    n=norm(x)
    if n == 0: raise ValueError('zero vector')
    return [complex(a/n) for a in x]

def eye(d): return [[complex(i==j) for j in range(d)] for i in range(d)]
def adj(A): return [list(map(complex.conjugate,row)) for row in zip(*A)]
def mm(A,B): return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]
def mv(A,x): return [sum(a*b for a,b in zip(row,x)) for row in A]
def sub(A,B): return [[a-b for a,b in zip(r,s)] for r,s in zip(A,B)]
def frob(A): return math.sqrt(sum(abs(z)**2 for row in A for z in row))
def outer(x,y): return [[a*b.conjugate() for b in y] for a in x]
def tensor(x,y): return [a*b for a in x for b in y]

def align(phi,chi):
    z=inner(phi,chi)
    phase=z.conjugate()/abs(z) if abs(z)>0 else 1
    return [phase*x for x in chi]

def chordal(phi,chi):
    return norm([a-b for a,b in zip(phi,align(phi,chi))])

def correction(phi,chi):
    """The explicit rank-two rotation; chi is phase aligned internally."""
    phi,chi=normalise(phi),normalise(chi)
    chi=align(phi,chi)
    if norm([a-b for a,b in zip(phi,chi)])<1e-13:
        return eye(len(phi)),chi
    c=max(0.0,min(1.0,inner(phi,chi).real))
    w=[b-c*a for a,b in zip(phi,chi)]
    s=norm(w)
    e=[z/s for z in w]
    pp,ee,ep,pe=outer(phi,phi),outer(e,e),outer(e,phi),outer(phi,e)
    I=eye(len(phi))
    R=[[I[i][j]+(c-1)*(pp[i][j]+ee[i][j])+s*(ep[i][j]-pe[i][j])
        for j in range(len(phi))] for i in range(len(phi))]
    return R,chi

def conditional(p):
    a=p[0]+p[1]
    return (a,)+((p[0]/a,) if a else ())+((p[2]/(1-a),) if a!=1 else ())

def admissible(p,L,occupation=False):
    return (sum(p)==1 and all(x>=0 for x in p)
            and all((L*x).denominator==1 for x in conditional(p))
            and (not occupation or all((L*x).denominator==1 for x in p)))

def probability_grid(A,B):
    return {(a*b,a*(1-b),(1-a)*c,(1-a)*(1-c))
        for m,n,r in product(range(A+1),range(B+1),range(B+1))
        for a,b,c in [(F(m,A),F(n,B),F(r,B))]}

def allocation(L,grid):
    assert L>=4 and L&(L-1)==0
    if grid=='T': return L,L
    M=L.bit_length()-1
    A=2**(M//2)
    return A,L//A

def rounded_prob(x,D): return min(D,max(0,math.floor(x*D+0.5)))

def phase_indices(phi,L):
    ref=next((z for z in phi if abs(z)>0),1)
    return [math.floor(cmath.phase(z/ref)*L/(2*math.pi)+0.5)%L if abs(z)>0 else 0
            for z in phi]

def reconstruct(probabilities,phases,L):
    return [math.sqrt(float(p))*cmath.exp(2j*math.pi*q/L)
            for p,q in zip(probabilities,phases)]

def quantise(phi,L,grid):
    """Return a numerical vector AND exact probability/phase certificates.

    This applies the deterministic coordinate quantiser. It does not attempt
    to decide exact membership of arbitrary floating-point input vectors.
    An exactly known admissible nominal output can instead be retained.
    """
    phi=normalise(phi)
    p=[abs(z)**2 for z in phi]
    total=sum(p);p=[x/total for x in p]
    A,B=allocation(L,grid)
    a=p[0]+p[1]; other=p[2]+p[3]
    b=p[0]/a if a>0 else 0
    c=p[2]/other if other>0 else 0
    m,n,r=rounded_prob(a,A),rounded_prob(b,B),rounded_prob(c,B)
    aa,bb,cc=F(m,A),F(n,B),F(r,B)
    probs=(aa*bb,aa*(1-bb),(1-aa)*cc,(1-aa)*(1-cc))
    phases=phase_indices(phi,L)
    assert admissible(probs,L,occupation=(grid=='K'))
    return reconstruct(probs,phases,L),{'A':A,'B':B,'m':m,'n':n,'r':r,
        'probabilities':list(map(str,probs)),'phase_indices':phases}

def quantise_local(phi,D,L):
    phi=normalise(phi)
    m=rounded_prob(abs(phi[0])**2,D)
    probs=(F(m,D),F(D-m,D))
    return reconstruct(probs,phase_indices(phi,L),L),probs

def eta(L,grid):
    A,B=allocation(L,grid)
    return math.sqrt(1/A+1/B)+2*math.sin(math.pi/(2*L))

def random_unitary(d,rng):
    U=eye(d)
    for i in range(d):
        for j in range(i+1,d):
            theta,alpha=rng.uniform(-math.pi,math.pi),rng.uniform(-math.pi,math.pi)
            c,s=math.cos(theta),math.sin(theta)
            G=eye(d)
            G[i][i]=G[j][j]=c
            G[i][j]=-cmath.exp(-1j*alpha)*s
            G[j][i]=cmath.exp(1j*alpha)*s
            U=mm(G,U)
    assert frob(sub(mm(adj(U),U),eye(d)))<TOL
    return U

def check_rotation(phi,chi):
    R,target=correction(phi,chi)
    unitary_error=frob(sub(mm(adj(R),R),eye(len(phi))))
    mapping_error=norm([a-b for a,b in zip(mv(R,phi),target)])
    # The two nonzero singular values of R-I are equal in the proof.
    norm_identity_error=abs(frob(sub(R,eye(len(phi))))/math.sqrt(2)-chordal(phi,chi))
    assert max(unitary_error,mapping_error,norm_identity_error)<TOL
    return unitary_error,mapping_error,norm_identity_error

def reduced_B(v):
    return [[sum(v[2*a+j]*v[2*a+k].conjugate() for a in range(2)) for k in range(2)] for j in range(2)]

def inv2(A):
    a,b=A[0];c,d=A[1];det=a*d-b*c
    return [[d/det,-b/det],[-c/det,a/det]]

def spectral_norm2(A):
    B=mm(adj(A),A)
    trace=(B[0][0]+B[1][1]).real
    discr=max(0.0,(B[0][0].real-B[1][1].real)**2+4*abs(B[0][1])**2)
    return math.sqrt(max(0.0,(trace+math.sqrt(discr))/2))

def verify():
    rng=random.Random(SEED)
    result={'seed':SEED,'numerical_tolerance':TOL,
        'exact_checks':{},'numerical_checks':{},
        'scope':'Exact certificates concern probabilities and integer phase indices. Matrix checks are floating point, not interval certificates.'}
    allocation_rows=[]
    for L in (4,8,16):
        T=probability_grid(L,L)
        K={p for p in T if admissible(p,L,True)}
        union=set()
        for r in range(L.bit_length()):
            A=2**r;union.update(probability_grid(A,L//A))
        assert union==K
        allocation_rows.append({'L':L,'K_probability_points':len(K),'allocation_union_points':len(union)})
    result['exact_checks']['allocation_decomposition']=allocation_rows
    budget_checks=0
    for L in (4,8,16,64):
        for m,n in product(range(L+1),repeat=2):
            a,b=F(m,L),F(n,L)
            p=(a*b,a*(1-b),(1-a)*b,(1-a)*(1-b))
            predicted=(L%(a.denominator*b.denominator)==0)
            assert admissible(p,L,True)==predicted
            budget_checks+=1
    result['exact_checks']['product_budget_cases']=budget_checks
    local_rows=[]
    for L in (4,8,16,64,256):
        b=F(L//2+1,L)
        allowed=[]
        for m in range(L+1):
            a=F(m,L);p=(a*b,a*(1-b),(1-a)*b,(1-a)*(1-b))
            if admissible(p,L,True):allowed.append(str(a))
        assert allowed==['0','1']
        local_rows.append({'L':L,'fixed_B_probability':str(b),'allowed_A_probabilities':allowed})
    result['exact_checks']['local_obstruction']=local_rows
    result['exact_checks']['local_operator_gap_formula']='sqrt(2-sqrt(2))'

    max_errors=[0.0,0.0,0.0]
    coverage=[]
    for L in (4,8,16,64,256,4096):
        for grid in ('T','K'):
            largest=0.0
            targets=[[complex(i==j) for i in range(4)] for j in range(4)]
            targets += [normalise([complex(rng.gauss(0,1),rng.gauss(0,1)) for _ in range(4)]) for _ in range(64)]
            for phi in targets:
                chi,cert=quantise(phi,L,grid)
                dist=chordal(phi,chi)
                assert dist<=eta(L,grid)+TOL
                largest=max(largest,dist)
                errs=check_rotation(phi,chi)
                max_errors=[max(a,b) for a,b in zip(max_errors,errs)]
            coverage.append({'L':L,'grid':grid,'targets_checked':len(targets),
                'largest_observed_chordal_distance':largest,'proved_uniform_bound':eta(L,grid)})
    check_rotation([1+0j,0j,0j,0j],[0j,1+0j,0j,0j])
    result['numerical_checks']['coverage']=coverage
    result['numerical_checks']['rotation_max_errors']=dict(zip(
        ('unitarity_frobenius','mapping_vector_norm','rotation_norm_identity'),max_errors))
    words=[]
    for L,grid in ((256,'T'),(4096,'K')):
        actual,_=quantise(normalise([1+2j,3-1j,-2j,1j]),L,grid)
        ideal=actual[:]
        cumulative=0.0
        for _ in range(24):
            U=random_unitary(4,rng)
            phi=mv(U,actual)
            target,_=quantise(phi,L,grid)
            R,aligned=correction(phi,target)
            Utilde=mm(R,U)
            assert chordal(mv(Utilde,actual),target)<TOL
            actual=aligned;ideal=mv(U,ideal);cumulative+=eta(L,grid)
            assert chordal(actual,ideal)<=cumulative+TOL
        words.append({'L':L,'grid':grid,'steps':24,'final_chordal_distance':chordal(actual,ideal),
            'accumulated_bound':cumulative})
    result['numerical_checks']['selected_words']=words

    # Product-sector construction: only the active factor is changed.
    L,A,B=256,16,16
    va,pa=quantise_local(normalise([1+1j,2-1j]),A,L)
    vb,pb=quantise_local(normalise([2-1j,-1+3j]),B,L)
    max_local_error=0.0
    for step in range(20):
        U=random_unitary(2,rng)
        if step%2==0:
            before=vb[:];phi=mv(U,va)
            va,pa=quantise_local(phi,A,L)
            check_rotation(phi,va)
            assert vb==before
            D=A
        else:
            before=va[:];phi=mv(U,vb)
            vb,pb=quantise_local(phi,B,L)
            check_rotation(phi,vb)
            assert va==before
            D=B
        selected=va if step%2==0 else vb
        err=chordal(phi,selected)
        assert err<=1/math.sqrt(D)+2*math.sin(math.pi/(2*L))+TOL
        max_local_error=max(max_local_error,err)
        assert all((L*x*y).denominator==1 for x in pa for y in pb)
    result['numerical_checks']['product_local_sector']={'L':L,'A':A,'B':B,'steps':20,
        'largest_observed_step_distance':max_local_error,'inactive_factor_unchanged':True}

    # Full-rank entangled fixed-marginal construction, coefficient matrices.
    Fmat=[[math.sqrt(3)/2,0j],[0j,0.5]]
    U=random_unitary(2,rng);Fmat=mm(U,Fmat)
    W=[[1+0j,0j],[0j,cmath.exp(0.17j)]]
    Cmat=mm(W,Fmat)
    phi=[z for row in Fmat for z in row];chi=[z for row in Cmat for z in row]
    aligned=align(phi,chi);Cmat=[aligned[:2],aligned[2:]]
    recovered=mm(Cmat,inv2(Fmat))
    assert frob(sub(reduced_B(phi),reduced_B(aligned)))<TOL
    assert frob(sub(mm(adj(recovered),recovered),eye(2)))<TOL
    bound=chordal(phi,aligned)/math.sqrt(0.25)
    observed=spectral_norm2(sub(recovered,eye(2)))
    assert observed<=bound+TOL
    result['numerical_checks']['fixed_marginal_formula']={'lambda_min':0.25,
        'observed_local_operator_distance':observed,'proved_candidate_bound':bound}
    result['status']='all checks passed'
    return result

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('nominal_gate_results.json'))
    args=parser.parse_args()
    result=verify()
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Exact allocation decomposition: K counts 27, 77, 233 at L=4,8,16.')
    print('Exact product denominator and local obstruction checks passed.')
    print('Numerical covering, unitary correction, circuit and local-sector checks passed.')
    print('Maximum rotation residuals:',result['numerical_checks']['rotation_max_errors'])
    print(result['status']);print('Wrote',args.output)

if __name__=='__main__':main()
