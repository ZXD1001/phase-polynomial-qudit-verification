#!/usr/bin/env python3
"""Independent native-Z6 probability reconstruction and local feedback audit.

No CRT or graph-stabilizer arithmetic is used to compute amplitudes. The exact
local inequalities are optimized over every second-layer decision for each
first-layer record. This is a finite regression, not the general proof.
"""
from __future__ import annotations
from itertools import product,combinations
from fractions import Fraction as F
from math import lcm
from pathlib import Path
import argparse,json,time
import numpy as np
if not __debug__:raise RuntimeError('Run without -O/-OO.')
ROOT=np.array([[1,0],[0,1],[-1,1],[-1,0],[0,-1],[1,-1]],dtype=np.int64)
BASE=Path(__file__).resolve().parents[1]

def exact_native(n):
    X=np.array(list(product(range(6),repeat=n)),dtype=np.int64)
    phase=((2*X[:,:-1]**2+3*X[:,:-1])*X[:,1:]).sum(axis=1)%6
    a0=ROOT[phase,0].reshape((6,)*n);b0=ROOT[phase,1].reshape((6,)*n)
    den=6**(2*n);dist={};power=6**np.arange(n-1,-1,-1)
    for setting in product(range(4),repeat=n):
        aa,bb=a0.copy(),b0.copy();nf=0
        for k,code in enumerate(setting):
            if not code:continue
            nf+=1;a=np.moveaxis(aa,k,0);b=np.moveaxis(bb,k,0);ca=np.zeros_like(a);cb=np.zeros_like(b)
            for j,x in product(range(6),repeat=2):
                c,d=ROOT[(-2*(code-1)*x*x+j*x)%6]
                ca[j]+=c*a[x]-d*b[x];cb[j]+=d*a[x]+(c+d)*b[x]
            aa=np.moveaxis(ca,0,k);bb=np.moveaxis(cb,0,k)
        N=(aa*aa+aa*bb+bb*bb).ravel();d0=6**(n+nf)
        assert int(N.sum())==d0
        dist[setting]=N*(den//d0)
    return X,power,den,dist

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--n',type=int,default=4)
    ap.add_argument('--output',type=Path,required=True);arg=ap.parse_args();n=arg.n
    if not 3<=n<=4:raise ValueError('This independent regression is limited to n=3,4.')
    start=time.monotonic();X,power,den,dist=exact_native(n);reject={}
    for s,p in dist.items():
        vv=[]
        for v in range(n):
            shift=np.zeros(n,dtype=int);shift[v]=2 if s[v] else 0
            index=((X+shift)%6)@power
            vv.append((p[index]*(p==0)).reshape((6,)*n))
        reject[s]=vv
    checks=0;type_counts={'late':0,'early_Z':0,'early_chirp_no_early_chirp_neighbor':0,'trivial_early_pair':0}
    maximum=F(-10);rows=[]
    for width in range(n):
      for A in combinations(range(n),width):
        B=tuple(j for j in range(n) if j not in A)
        for front in product(range(4),repeat=len(B)):
            fixed=dict(zip(B,front));best=[np.full(6**len(B),-den*10,dtype=np.int64) for _ in range(n)]
            kinds=[]
            for v in range(n):
                neighbors=[j for j in (v-1,v+1) if 0<=j<n]
                if v in A:k='late'
                elif fixed[v]==0:k='early_Z'
                elif any(j in B and fixed[j]!=0 for j in neighbors):k='trivial_early_pair'
                else:k='early_chirp_no_early_chirp_neighbor'
                kinds.append(k)
            for terminal in product(range(4),repeat=width):
                s=[0]*n
                for j,c in zip(B,front):s[j]=c
                for j,c in zip(A,terminal):s[j]=c
                s=tuple(s);pf=dist[s].reshape((6,)*n).sum(axis=A).ravel()
                for v in range(n):
                    q=reject[s][v].sum(axis=A).ravel()
                    if kinds[v]=='late':z=q-int(s[v]!=0)*pf
                    elif kinds[v]=='early_Z':z=q
                    elif kinds[v]=='trivial_early_pair':z=q
                    else:
                        bad=any(s[j]!=0 for j in (v-1,v+1) if 0<=j<n)
                        z=q-int(bad)*pf
                    np.maximum(best[v],z,out=best[v])
            for v,k in enumerate(kinds):
                # The final vertex has no quadratic coefficient and is not part of
                # the new qutrit-error sum. We only assert its late/Z facts here.
                if v==n-1 and k.startswith('early_chirp'):continue
                val=F(int(best[v].sum()),den)
                limit=F(0) if k in ('late','early_Z') else F(1) if k=='trivial_early_pair' else F(1,3)
                assert val<=limit,(n,A,front,v,k,val,limit)
                if k=='early_chirp_no_early_chirp_neighbor':maximum=max(maximum,val)
                type_counts[k]+=1;checks+=1
    # Direct independent check of the full positive witness, at width one.
    cert=json.loads((BASE/'inputs'/f'robust_n{n}_w1.json').read_text())
    weights=[F(t['weight']) for t in cert['terms']];wd=lcm(*(w.denominator for w in weights));wi=[int(w*wd) for w in weights]
    errors=np.array([(3*np.array(t['qubit_character'])+2*np.array(t['qutrit_character']))%6 for t in cert['terms']])
    maxarr=np.zeros((n*4**(n-1),6**(n-1)),dtype=np.int64)
    for s,p in dist.items():
        idx=((X[None,:,:]+errors[:,None,:]*np.array(s,dtype=bool)[None,None,:])%6)@power
        q=(np.array(wi)[:,None]*p[idx]).sum(axis=0)*(p==0)
        assert int(q.sum())<=den*wd;q=q.reshape((6,)*n)
        for k in range(n):
            front=s[:k]+s[k+1:];j=k*4**(n-1)+sum(c*4**(n-2-i) for i,c in enumerate(front))
            np.maximum(maxarr[j],q.sum(axis=k).ravel(),out=maxarr[j])
    totals=maxarr.sum(axis=1);mx=F(int(totals.max()),den*wd)
    old=json.loads((BASE/'results'/f'robust_n{n}_w1_native_verified.json').read_text())
    assert all(F(int(x),den*wd)==F(y,old['common_denominator']) for x,y in zip(totals,old['per_template_maximum_numerators']))
    out=dict(status='EXACT_PASS',n=n,seconds=time.monotonic()-start,
        local_constraints=checks,local_constraint_types=type_counts,
        maximum_repair_penalized_rejection=str(maximum),repair_upper='1/3',
        witness_maximum=str(mx),all_template_values_equal_independent_CRT=True,
        native_templates=len(totals),arithmetic='Z[omega6] integer pairs; exact squared norms, zero support, sums and per-record maxima.',
        note='Optimizes every late basis choice for each raw six-valued front record; random mixtures are covered by convexity.')
    arg.output.parent.mkdir(parents=True,exist_ok=True);arg.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
