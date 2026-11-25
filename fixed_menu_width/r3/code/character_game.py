"""Cutting-plane minimization over general CRT character-error mixtures.
Exact support and fixed-basis acceptance are constructed from integer data.
Floating LP proposes witnesses; a separate verifier certifies them.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
from itertools import product
import numpy as np,json,time,sys
from scipy.optimize import linprog
from menu_witness import configs,distributions
BASE=Path(__file__).resolve().parents[1]

def qubit_accept_table(n):
    X=configs(2,n);adj=np.c_[np.zeros(2**n,dtype=int),X[:,:-1]]^np.c_[X[:,1:],np.zeros(2**n,dtype=int)];out={}
    for axes in product(range(2),repeat=n):
        S=X[np.all(np.where(np.array(axes)[None,:],adj,X)==0,axis=1)]
        out[axes]=np.all((X@S.T)%2==0,axis=1).astype(int)
    return X,out

def error_set(n,mode):
    X2=configs(2,n);X3=configs(3,n)
    if mode=='full':return [(tuple(a),tuple(b)) for a in X2 for b in X3 if any(a) or any(b)]
    z=(0,)*n
    one3=[]
    for v in range(n):
        for t in [1,2]:
            e=[0]*n;e[v]=t;one3.append(tuple(e))
    if mode=='local':return [(z,b) for b in one3]+[(tuple(a),z) for a in X2[1:]]
    if mode=='localmixed':return [(tuple(a),b) for a in X2 for b in [z]+one3 if any(a) or any(b)]
    raise ValueError(mode)

def build(n,mode):
    tic=time.time();X3,probs,nums,dens=distributions(n);X2,a2=qubit_accept_table(n);errs=error_set(n,mode)
    e2=np.array([a for a,b in errs]);e3=np.array([b for a,b in errs]);p2=2**np.arange(n-1,-1,-1);p3=3**np.arange(n-1,-1,-1)
    i2=e2@p2;ng=n*4**(n-1);no=3**(n-1);ne=len(errs)
    G=np.empty((ng,no,4,ne),dtype=float)
    for s,prob in probs.items():
        ax=np.array(s,dtype=bool)
        i3=((X3[None,:,:]-e3[:,None,:]*ax[None,None,:])%3)@p3
        r=prob[i3]*(1-a2[tuple(ax.astype(int))][i2,None]*(nums[s]!=0)[None,:])
        for k in range(n):
            front=s[:k]+s[k+1:];idx=k*4**(n-1)+int(np.dot(front,4**np.arange(n-2,-1,-1)))
            G[idx,:,s[k],:]=r.reshape((ne,)+(3,)*n).sum(axis=k+1).reshape(ne,-1).T
    labels=[(k,s) for k in range(n) for s in product(range(4),repeat=n-1)]
    print('BUILT',n,mode,G.shape,'seconds',time.time()-tic,flush=True)
    return G,errs,labels

def solve(n,mode='full',maxit=200):
    G,errs,labels=build(n,mode);ng,no,nb,ne=G.shape
    cuts=[];w=np.ones(ne)/ne;best=1.;bestw=w.copy();hist=[];tic=time.time()
    for it in range(maxit):
        rr=np.einsum('gobe,e->gob',G,w,optimize=True)
        choice=rr.argmax(axis=2);score=rr.max(axis=2).sum(axis=1);worst=float(score.max())
        if worst<best:best=worst;bestw=w.copy()
        # Add several most violated trees; all are genuine members of the class.
        ids=np.argsort(score)[-min(32,ng):]
        for g in ids:cuts.append(G[g,np.arange(no),choice[g],:].sum(axis=0))
        C=np.array(cuts)
        lp=linprog(np.r_[np.zeros(ne),1.],A_ub=np.c_[C,-np.ones(len(C))],b_ub=np.zeros(len(C)),A_eq=[np.r_[np.ones(ne),0.]],b_eq=[1.],bounds=[(0,None)]*(ne+1),method='highs')
        if not lp.success:raise RuntimeError(lp.message)
        w=lp.x[:ne];lo=float(lp.x[-1]);hist.append([it,lo,best])
        print('iter',it,'diagonal_game_lower',lo,'valid_witness_upper',best,'time',time.time()-tic,flush=True)
        if best-lo<2e-9:break
        if len(cuts)>1000:
            active=np.flatnonzero(-lp.ineqlin.marginals>1e-10).tolist();keep=sorted(set(active+list(range(len(cuts)-200,len(cuts)))))
            cuts=[cuts[j] for j in keep]
    from fractions import Fraction
    support=[{'weight_float':float(v),'weight_rational_proposal':str(Fraction(float(v)).limit_denominator(1000000)),'qubit_character':list(map(int,a)),'qutrit_character':list(map(int,b))} for (a,b),v in zip(errs,bestw) if v>1e-9]
    out={'n':n,'mode':mode,'status':'FLOATING_CHARACTER_GAME','witness_upper':best,'diagonal_game_lower_NOT_quantum_lower':lo,'support':support,'history':hist}
    (BASE/'results'/f'character_game_n{n}_{mode}.json').write_text(json.dumps(out,indent=2))
    print('SUPPORT',json.dumps(support),flush=True)
    return out
if __name__=='__main__':solve(int(sys.argv[1]) if len(sys.argv)>1 else 4,sys.argv[2] if len(sys.argv)>2 else 'full')
