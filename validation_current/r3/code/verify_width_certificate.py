"""Verify finite same-menu adaptive-width upper bounds using exact arithmetic.

--native: exhaustive raw 2^n x 3^n product distributions (no twirl/reduction).
Default: exact CRT character-error reduction proved in the accompanying report.
Each invocation enumerates subsets of EXACTLY the indicated size. To bound
at-most r sites, combine checks of sizes 1,...,r. Size zero embeds in size one.
All second-layer basis choices are optimized separately at each first-layer record,
so arbitrary correlations and record dependence are included in the upper bound.
"""
if not __debug__:
    raise RuntimeError("Exact verification requires Python without -O/-OO.")

import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
from fractions import Fraction as F
from itertools import product,combinations
from math import lcm,comb
import numpy as np,json,time,sys,hashlib
from menu_witness import distributions,configs
from character_game import qubit_accept_table
from verify_native_certificate import qubit_dist
BASE=Path(__file__).resolve().parents[1]

def verify(path,native=False,strict=True):
    tic=time.time();path=Path(path);cert=json.loads(path.read_text());n=cert['n'];w=cert['adaptive_width'];terms=cert['terms'];claim=F(cert['claim_upper'])
    assert isinstance(n,int) and isinstance(w,int) and 2<=n and 0<=w<n
    assert cert['menu']=='Z,B0,B1,B2' and terms
    weights=[F(t['weight']) for t in terms];assert sum(weights)==1 and all(t>=0 for t in weights)
    wd=lcm(*(v.denominator for v in weights));wi=np.array([int(v*wd) for v in weights],dtype=np.int64)
    e2=np.array([t['qubit_character'] for t in terms],int);e3=np.array([t['qutrit_character'] for t in terms],int)
    assert all(any(a) or any(b) for a,b in zip(e2,e3))
    assert e2.shape==e3.shape==(len(terms),n) and np.all((e2>=0)&(e2<2)) and np.all((e3>=0)&(e3<3))
    X3,probs,N3,D3=distributions(n);del probs
    X2,a2=qubit_accept_table(n)
    if native:X2,N2=qubit_dist(n)
    p2=2**np.arange(n-1,-1,-1);p3=3**np.arange(n-1,-1,-1)
    subs=list(combinations(range(n),w));ng=len(subs)*4**(n-w);nr=(6 if native else 3)**(n-w)
    maxima=np.zeros((ng,nr),dtype=np.int64)
    den=wd*((6 if native else 3)**(2*n))
    assert den < np.iinfo(np.int64).max//4,'Use Python integer arrays for larger inputs.'
    for s,tar3 in N3.items():
        ax=np.array(s,dtype=bool);ix3=((X3[None,:,:]-e3[:,None,:]*ax[None,None,:])%3)@p3
        if native:
            tar2,d2=N2[tuple(ax.astype(int))];r=np.zeros((2**n,3**n),dtype=np.int64)
            ix2=((X2[None,:,:]-e2[:,None,:]*ax[None,None,:])%2)@p2
            for j,ww in enumerate(wi):r+=ww*tar2[ix2[j],None]*tar3[None,ix3[j]]
            assert int(r.sum())==wd*d2*D3[s]
            r*=~((tar2[:,None]!=0)&(tar3[None,:]!=0));r*=6**(2*n)//(d2*D3[s]);r=r.reshape((2,)*n+(3,)*n)
        else:
            acc=a2[tuple(ax.astype(int))][e2@p2]
            r=(wi[:,None]*tar3[ix3]*(1-acc[:,None]*(tar3!=0)[None,:])).sum(axis=0)
            r*=3**(2*n)//D3[s];r=r.reshape((3,)*n)
        assert int(r.sum())<=den
        for ia,A in enumerate(subs):
            front=tuple(s[j] for j in range(n) if j not in A)
            ig=ia*4**(n-w)+int(np.dot(front,4**np.arange(n-w-1,-1,-1)))
            axes=A+tuple(n+j for j in A) if native else A
            v=r.sum(axis=axes).ravel()
            np.maximum(maxima[ig],v,out=maxima[ig])
    sums=maxima.sum(axis=1);mx=int(sums.max());up=F(mx,den);passed=up<=claim
    bysubset=[]
    for ia,A in enumerate(subs):
        q=F(int(sums[ia*4**(n-w):(ia+1)*4**(n-w)].max()),den)
        bysubset.append({'second_layer_sites':A,'maximum_rejection':str(q)})
    result={'status':'EXACT_PASS' if passed else 'CLAIM_FAILED_EXACT_EVALUATION','n':n,'adaptive_width':w,'claim_upper':str(claim),
      'actual_maximum_rejection':str(up),'decimal_display':float(up),'native_exhaustive':native,'enumerated_exact_second_layer_size':w,'at_most_width_note':'Combine every positive size up to requested width; nonadaptive size zero embeds in size one.','templates':ng,'front_records_per_template':nr,
      'second_layer_choices_per_record':4**w,'full_basis_patterns':4**n,'witness_terms':len(terms),'common_denominator':den,'seconds':time.time()-tic,
      'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'per_subset':bysubset,
      'per_template_maximum_numerators':sums.tolist(),'arithmetic':'Integer cyclotomic transforms and squared norms; exact support; max and sum over all branches; no eigenvalue or floating feasibility check.'}
    dest=BASE/'results'/(path.stem+('_native' if native else '_reduced')+'_verified.json');dest.write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k not in ['per_subset','per_template_maximum_numerators']}),flush=True)
    if strict and not passed:raise AssertionError(f'Claim {claim} FAILED; actual {up}.')
    return result
if __name__=='__main__':
    for arg in sys.argv[1:]:
        if not arg.startswith('--'):verify(arg,native='--native' in sys.argv,strict='--allow-failure' not in sys.argv)
