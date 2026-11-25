#!/usr/bin/env python3
"""Exact finite checks for the path X/Z support and strict dual-witness lemma.

Independent GF(2) row reduction versus the odd-X-run rule. Checks all X/Z
patterns for n=2,...,12 by default. Optionally tests the leaf-overlap identity
by direct integer Walsh transforms for all outcomes/errors through n=6.
These tests are not an infinite-size proof; REPORT.md provides that proof.
"""
from __future__ import annotations

if not __debug__:
    raise RuntimeError("Run exact checks without -O or -OO.")
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import argparse, json, time
import numpy as np


def nullspace(rows: list[int],n:int)->list[int]:
    rr=list(rows);pivots=[];i=0
    for col in range(n):
        j=next((j for j in range(i,len(rr)) if (rr[j]>>col)&1),None)
        if j is None:continue
        rr[i],rr[j]=rr[j],rr[i]
        for j in range(len(rr)):
            if j!=i and (rr[j]>>col)&1:rr[j]^=rr[i]
        pivots.append(col);i+=1
    basis=[]
    for free in set(range(n))-set(pivots):
        v=1<<free
        for row,col in zip(rr[:i],pivots):
            if (row>>free)&1:v|=1<<col
        basis.append(v)
    return basis


def run_basis(axes:tuple[int,...])->list[int]:
    n=len(axes);out=[];v=0
    while v<n:
        if not axes[v]:v+=1;continue
        end=v
        while end<n and axes[end]:end+=1
        if (end-v)%2:out.append(sum(1<<j for j in range(v,end,2)))
        v=end
    return out


def span(basis:list[int])->set[int]:
    out={0}
    for a in basis:out|={v^a for v in list(out)}
    return out


def reject(basis:list[int],error:int)->int:
    return int(any((a&error).bit_count()%2 for a in basis))


def check_patterns(max_n:int)->dict:
    start=time.monotonic();rows_out=[]
    for n in range(2,max_n+1):
        den=4*n-6;max_noncolor=F(0);n_equal=0;checked=0
        for axes in product((0,1),repeat=n):
            rows=[((1<<(j-1) if j else 0)|(1<<(j+1) if j+1<n else 0))
                  if axes[j] else (1<<j) for j in range(n)]
            a=nullspace(rows,n);b=run_basis(axes)
            assert span(a)==span(b)
            I=0
            for v in b:I|=v
            assert (I & (I<<1))==0
            numerator=0
            for v in range(n-1):
                z=reject(b,1<<v)+reject(b,1<<(v+1))
                assert z<=1
                numerator+=z
            for v in range(n-2):
                z=reject(b,(1<<v)|(1<<(v+2)))+reject(b,1<<(v+1))
                assert z<=1
                numerator+=z
            is_color=all(axes[j]!=axes[j+1] for j in range(n-1))
            assert (numerator==2*n-3)==is_color
            if not is_color:
                assert numerator<=2*n-4
                max_noncolor=max(max_noncolor,F(numerator,den))
            else:n_equal+=1
            checked+=1
        assert n_equal==2
        rows_out.append({'n':n,'patterns':checked,'half_equality_patterns':n_equal,
                         'max_noncolor_rejection':str(max_noncolor),
                         'proved_general_upper_noncolor':str(F(1,2)-F(1,den))})
    return {'status':'EXACT_PASS','check':'Independent GF2 elimination, run decomposition, '
        'local edge/triple pair constraints, and iff color equality',
        'seconds':time.monotonic()-start,'cases':rows_out}


def check_walsh(max_n:int)->dict:
    start=time.monotonic();out=[]
    for n in range(2,max_n+1):
        X=np.array(list(product((0,1),repeat=n)),dtype=np.int64)
        base=(1-2*((X[:,:-1]*X[:,1:]).sum(axis=1)%2)).reshape((2,)*n)
        pow2=2**np.arange(n-1,-1,-1);pairs=0
        for axes in product((0,1),repeat=n):
            aa=base.copy()
            for j,ax in enumerate(axes):
                if ax:
                    vv=np.moveaxis(aa,j,0)
                    aa=np.moveaxis(np.stack((vv[0]+vv[1],vv[0]-vv[1])),0,j)
            N=(aa*aa).ravel();active=np.flatnonzero(N);k=len(run_basis(axes))
            assert len(active)==2**(n-k)
            assert np.all(N[active]==N[active[0]])
            for e in X:
                mask=sum(int(e[j])<<j for j in range(n))
                index=((X-e*np.array(axes))%2)@pow2
                expected=N[active]*int(not reject(run_basis(axes),mask))
                assert np.array_equal(N[index[active]],expected)
                pairs+=len(active)
        out.append({'n':n,'basis_patterns':2**n,'positive_leaf_error_comparisons':pairs})
    return {'status':'EXACT_PASS','check':'Direct integer Walsh positive-leaf likelihood-ratio identity',
            'seconds':time.monotonic()-start,'cases':out}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--max-n',type=int,default=12)
    p.add_argument('--walsh-max-n',type=int,default=6)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();out={'patterns':check_patterns(a.max_n),'walsh':check_walsh(a.walsh_max_n)}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
