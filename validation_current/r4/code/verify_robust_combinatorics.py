#!/usr/bin/env python3
"""Finite exact tests for the new counting lemmas; NOT an all-n proof.

Uses independent GF(2) elimination for S_tau, rather than the run formula.
Tests full X/Z types through max_n, all early/late masks through mask_max_n,
and derives near-optimal integer resource tables by exact Fraction arithmetic.
"""
from __future__ import annotations
import argparse,json,time
from itertools import product
from pathlib import Path
from fractions import Fraction as F
from verify_graph_lemmas import nullspace,reject,run_basis,span
from robust_certificate import generate
if not __debug__:raise RuntimeError('Run without -O/-OO.')

def data(tau):
    n=len(tau)
    rows=[((1<<(j-1) if j else 0)|(1<<(j+1) if j+1<n else 0)) if tau[j] else 1<<j for j in range(n)]
    basis=nullspace(rows,n);assert span(basis)==span(run_basis(tau))
    I=[reject(basis,1<<j) for j in range(n)]
    H=[int(tau[j] and ((j and tau[j-1]) or (j+1<n and tau[j+1]))) for j in range(n)]
    K=sum(tau)-sum(I)
    ds=[1-reject(basis,(1<<j)|(1<<(j+2)))-I[j+1] for j in range(n-2)]
    assert all(d in (0,1) for d in ds)
    return I,H,K,sum(ds)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--max-n',type=int,default=13);ap.add_argument('--mask-max-n',type=int,default=8)
    ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();start=time.monotonic();cases=[]
    for n in range(3,a.max_n+1):
        count=0;masks=0;tight=0
        for tau in product((0,1),repeat=n):
            I,H,K,T=data(tau)
            assert K+2*sum(H)<=7*T,(n,tau,I,H,K,T)
            tight+=int(K+2*sum(H)==7*T)
            if n%2:assert sum(I[:-1])<=F(n-1,2)
            else:assert sum(I[:-1])+I[-1]<=F(n,2)
            if n<=a.mask_max_n:
                for amask in range(1<<n):
                    early=[int(tau[j] and not(amask>>j&1)) for j in range(n)]
                    fixed=sum(1 for j in range(n) if early[j] and
                              ((j and early[j-1]) or (j+1<n and early[j+1])))
                    dynamic=sum(1 for j in range(n) if early[j] and H[j] and
                                not((j and early[j-1]) or (j+1<n and early[j+1])))
                    assert 2*fixed+3*dynamic<=2*sum(H),(n,tau,amask)
                    # Three times the pointwise correction term.
                    assert K+2*fixed+3*dynamic<=7*T
                    masks+=1
            count+=1
        cases.append(dict(n=n,patterns=count,pattern_mask_pairs=masks,count_bound_equalities=tight))
    table=[]
    for n in [10,100,1000,10000]:
        D=17*n-31+int(n%2==0)
        for Delta in [F(1,100),F(1,50),F(1,25),F(1,20)]:
            raw=F(n-1,2)-F(D,2)*Delta
            lower=max(0,-(-raw.numerator//raw.denominator))
            required=F(1,2)-Delta
            def feasible(r):
                m=[len(range(c,n-1,2)) for c in (0,1)]
                q=[F(1,3)+F(2*min(r,z),3*z) for z in m]
                return q[0]*q[1]/sum(q)
            upper=next(r for r in range(n//2+1) if feasible(r)>=required)
            assert lower<=upper
            table.append(dict(n=n,Delta=str(Delta),target_gap=str(required),
                              width_necessary=lower,width_sufficient=upper))
    result=dict(status='EXACT_PASS',seconds=time.monotonic()-start,cases=cases,
        inequalities=['K+2H <= 7T for all n>=3 (finite regression)',
                      '2F+3D_dynamic <= 2H for every early/late mask (finite regression)',
                      'parity-specific independent-set packing'],resource_table=table,
        note='Infinite-size claims follow from analytic proofs, not these finite enumerations.')
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
