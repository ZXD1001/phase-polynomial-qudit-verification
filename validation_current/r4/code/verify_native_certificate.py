"""Exact all-leaf verification of an adaptive same-menu witness, n=4.

All 4^n full bases, all 6^n outcomes, all one-terminal decision rules are
covered by summing a per-front-outcome maximum over four terminal bases.
Only bounded int64 arithmetic (bounds checked), with Fraction at output.
No optimizer, floating tolerance, or matrix eigenvalue is used.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from itertools import product
from fractions import Fraction
from pathlib import Path
import numpy as np,json,time
from menu_witness import distributions,configs
BASE=Path(__file__).resolve().parents[1]

def qubit_dist(n):
    X=configs(2,n);p=(X[:,:-1]*X[:,1:]).sum(axis=1)%2
    A=(1-2*p).reshape((2,)*n);out={}
    for axes in product(range(2),repeat=n):
        aa=A.copy();nf=sum(axes)
        for k,t in enumerate(axes):
            if not t:continue
            a=np.moveaxis(aa,k,0);b=np.stack([a[0]+a[1],a[0]-a[1]],axis=0);aa=np.moveaxis(b,0,k)
        N=(aa*aa).ravel();den=2**(n+nf)
        assert int(N.sum())==den
        out[axes]=(N,den)
    return X,out

def verify(n=4):
    tic=time.time();X3,_,N3,D3=distributions(n);X2,N2=qubit_dist(n)
    # Symmetry-averaged certificate: four qutrit errors weight 9/70 each.
    terms=[]
    for v in [0,2]:
        for t in [1,2]:
            e=[0]*n;e[v]=t;terms.append((9,[0]*n,e))
    for w,positions in [(5,[2]),(18,[1]),(11,[1,3])]:
        e=[0]*n
        for j in positions:e[j]=1
        terms.append((w,e,[0]*n))
    assert sum(w for w,_,_ in terms)==70
    den=70*6**(2*n);power2=2**np.arange(n-1,-1,-1);power3=3**np.arange(n-1,-1,-1)
    full={};maxnum=0
    for s in product(range(4),repeat=n):
        axes=tuple(int(t!=0) for t in s);a2,d2=N2[axes];a3=N3[s];d3=D3[s]
        target_support=(a2[:,None]!=0)&(a3[None,:]!=0)
        rej=np.zeros((2**n,3**n),dtype=np.int64)
        for w,e2,e3 in terms:
            i2=((X2-np.array(e2)*np.array(axes))%2)@power2
            i3=((X3-np.array(e3)*np.array(axes))%3)@power3
            scale=6**(2*n)//(d2*d3)
            rej+=w*scale*a2[i2,None]*a3[None,i3]
        assert int(rej.sum())==den
        rej*=~target_support
        full[s]=rej.reshape((2,)*n+(3,)*n)
        maxnum=max(maxnum,int(rej.max()))
    values=[]
    for k in range(n):
        for front in product(range(4),repeat=n-1):
            choices=[]
            for b in range(4):
                s=list(front);s.insert(k,b)
                # Sum over both CRT outcomes at the terminal.
                z=full[tuple(s)].sum(axis=n+k).sum(axis=k)
                choices.append(z)
            num=int(np.max(np.stack(choices),axis=0).sum())
            q=Fraction(num,den)
            values.append({'terminal':k,'front_string':front,'maximum_rejection':str(q),'maximizing_choices':np.argmax(np.stack(choices),axis=0).ravel().tolist()})
    val=max(Fraction(v['maximum_rejection']) for v in values)
    if n==4:assert val==Fraction(29,70),(val,float(val))
    result={'status':'EXACT_PASS','n':n,'class':'one prechosen terminal; all front and terminal PVMs in {Z,B0,B1,B2}; arbitrary front/terminal probabilities and arbitrary dependence on the full native-six-level front record; perfect completeness; single pass per site',
      'witness_terms':[{'weight':str(Fraction(w,70)),'qubit_character':e2,'qutrit_character':e3} for w,e2,e3 in terms],
      'upper_bound':str(val),'upper_bound_decimal_display':float(val),'total_front_templates':len(values),'native_front_records_per_template':6**(n-1),'terminal_bases_per_record':4,
      'common_probability_denominator':den,'maximum_leaf_numerator':maxnum,'seconds':time.time()-tic,'per_template':values,
      'PSD_orthogonality_reason':'Positive rational mixture of normalized nonzero character errors on constant-modulus pure factors; finite character orthogonality is exact.',
      'arithmetic':'integer cyclotomic amplitude transforms, exact squared norms, exact supports, exact max and sum. NumPy int64 entries bounded by common denominator; Fraction output.'}
    (BASE/'results'/f'exact_native_witness_n{n}.json').write_text(json.dumps(result,indent=2))
    print('EXACT',n,str(val),'templates',len(values),'seconds',time.time()-tic,flush=True)
    return result
if __name__=='__main__':
    import sys
    for n in map(int,sys.argv[1:] or [4]):verify(n)
