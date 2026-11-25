"""Independent n=4 NEW FAMILY certificate reconstruction directly in Z[omega_6].
No CRT amplitude factorization and no graph-stabilizer code is used.
The witness characters j=3a+2b mod6 are applied as native Z6 errors.
"""
if not __debug__:
    raise RuntimeError('Run without -O/-OO.')
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from math import lcm
import numpy as np,json,time
ROOT=np.array([[1,0],[0,1],[-1,1],[-1,0],[0,-1],[1,-1]],dtype=np.int64)
BASE=Path(__file__).resolve().parents[1]

def main():
    tic=time.time();n=4;X=np.array(list(product(range(6),repeat=n)),dtype=np.int64)
    phase=((2*X[:,:-1]**2+3*X[:,:-1])*X[:,1:]).sum(axis=1)%6
    a0=ROOT[phase,0].reshape((6,)*n);b0=ROOT[phase,1].reshape((6,)*n)
    cert=json.loads((BASE/'inputs'/'family_n4_w1.json').read_text())
    wt=[F(z['weight']) for z in cert['terms']];wd=lcm(*(w.denominator for w in wt));wi=[int(w*wd) for w in wt]
    errors=np.array([(3*np.array(z['qubit_character'])+2*np.array(z['qutrit_character']))%6 for z in cert['terms']])
    den=wd*6**(2*n);pow6=6**np.arange(n-1,-1,-1);maxima=np.zeros((n*4**(n-1),6**(n-1)),dtype=np.int64)
    for setting in product(range(4),repeat=n):
        aa,bb=a0.copy(),b0.copy();nf=0
        for k,code in enumerate(setting):
            if not code:continue
            nf+=1;a=np.moveaxis(aa,k,0);b=np.moveaxis(bb,k,0);ca=np.zeros_like(a);cb=np.zeros_like(b)
            for j,x in product(range(6),repeat=2):
                c,d=ROOT[(-2*(code-1)*x*x+j*x)%6]
                ca[j]+=c*a[x]-d*b[x];cb[j]+=d*a[x]+(c+d)*b[x]
            aa=np.moveaxis(ca,0,k);bb=np.moveaxis(cb,0,k)
        N=(aa*aa+aa*bb+bb*bb).ravel();d0=6**(n+nf);assert int(N.sum())==d0
        # Native inverse Fourier convention: Z6 error shifts the bra outcome by +e.
        ix=((X[None,:,:]+errors[:,None,:]*np.array(setting,dtype=bool)[None,None,:])%6)@pow6
        prob=(np.array(wi)[:,None]*N[ix]).sum(axis=0);assert int(prob.sum())==wd*d0
        prob*=N==0;prob*=6**(2*n)//d0;prob=prob.reshape((6,)*n)
        for k in range(n):
            front=setting[:k]+setting[k+1:];idx=k*4**(n-1)+sum(t*4**(n-2-j) for j,t in enumerate(front))
            np.maximum(maxima[idx],prob.sum(axis=k).ravel(),out=maxima[idx])
    totals=maxima.sum(axis=1);upper=F(int(totals.max()),den);assert upper==F(cert['claim_upper'])
    prior=json.loads((BASE/'results'/'family_n4_w1_native_verified.json').read_text())
    assert all(F(int(x),den)==F(y,prior['common_denominator']) for x,y in zip(totals,prior['per_template_maximum_numerators']))
    result={'status':'EXACT_PASS','n':4,'arithmetic':'Independent native-six-level integer cyclotomic amplitudes; omega6^2=omega6-1; no CRT factorization in amplitude computation',
      'upper':str(upper),'all_256_template_values_agree_with_CRT_verifier':True,'common_denominator':den,'seconds':time.time()-tic,
      'per_template_maximum_numerators':totals.tolist()}
    (BASE/'results'/'family_direct_omega6_n4_verified.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k!='per_template_maximum_numerators'}),flush=True)
if __name__=='__main__':main()
