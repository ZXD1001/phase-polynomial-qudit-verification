"""Same-menu one-terminal adaptive support tests: character-diagonal dual.

Qutrit front and terminal may use all four Pauli bases C0=Z,
C(1+t)[x,j] = omega**(t*x*x+j*x)/sqrt(3).
The dual is a mixture of exact orthogonal character errors Z^e |H_n>.
Each front configuration may choose its terminal basis on every record.
Zero support is decided by integer cyclotomic arithmetic, never tolerance.
LP discovery is floating point; rational certificates are verified separately.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
from itertools import product
import json,time,sys
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix,hstack,csr_matrix
ROOT=np.array([[1,0],[0,1],[-1,-1]],dtype=np.int64)
BASE=Path(__file__).resolve().parents[1]

def configs(d,n):return np.array(list(product(range(d),repeat=n)),dtype=np.int64)

def transform_pair(A,B,axis,t):
    # Bra of chirp basis; normalization omitted.
    aa=np.moveaxis(A,axis,0);bb=np.moveaxis(B,axis,0)
    oa=np.zeros_like(aa);ob=np.zeros_like(bb)
    for j in range(3):
        for x in range(3):
            c,d=ROOT[(-t*x*x-j*x)%3]
            oa[j]+=c*aa[x]-d*bb[x]
            ob[j]+=d*aa[x]+(c-d)*bb[x]
    return np.moveaxis(oa,0,axis),np.moveaxis(ob,0,axis)

def distributions(n):
    X=configs(3,n);ph=np.sum(X[:,:-1]**2*X[:,1:],axis=1)%3
    A=ROOT[ph,0].reshape((3,)*n);B=ROOT[ph,1].reshape((3,)*n)
    probs={};nums={};dens={}
    for s in product(range(4),repeat=n):
        aa,bb=A.copy(),B.copy();nf=0
        for k,t in enumerate(s):
            if t:aa,bb=transform_pair(aa,bb,k,t-1);nf+=1
        N=(aa*aa-aa*bb+bb*bb).ravel();den=3**(n+nf)
        assert int(N.sum())==den
        probs[s]=N/den;nums[s]=N;dens[s]=den
    return X,probs,nums,dens

def build(n,errors=None):
    tic=time.time();X,probs,nums,dens=distributions(n);D=3**n
    if errors is None:errors=X[1:]
    E=len(errors);powers=3**np.arange(n-1,-1,-1)
    # For each total basis s: R_s[e,o_full] exact support zero times error prob.
    rs={}
    for s in product(range(4),repeat=n):
        shifts=errors*np.array(s,dtype=bool)[None,:]
        indices=((X[None,:,:]-shifts[:,None,:])%3)@powers
        rs[s]=probs[s][indices]*(nums[s]==0)[None,:]
    groups=[];labels=[]
    for k in range(n):
        sites=[j for j in range(n) if j!=k]
        for front in product(range(4),repeat=n-1):
            R=[]
            for b in range(4):
                s=list(front);s.insert(k,b);arr=rs[tuple(s)].reshape((E,)+(3,)*n)
                R.append(arr.sum(axis=k+1).reshape(E,-1).T)
            groups.append(np.stack(R,axis=1));labels.append((k,front))
    # G[group,front_outcome,terminal_basis,error]
    G=np.stack(groups);print('build',n,G.shape,'seconds',time.time()-tic,flush=True)
    return G,labels,errors,nums,dens

def search(n, errors=None):
    G,labels,errors,nums,dens=build(n,errors);ng,no,nb,ne=G.shape
    ynum=ng*no
    # min t s.t. G[g,o,b,:]w <= y[g,o], sum_o y[g,o]<=t.
    row=np.arange(ng*no*nb);cols=np.repeat(np.arange(ynum),nb)
    Y=coo_matrix((-np.ones(len(row)),(row,cols)),shape=(len(row),ynum)).tocsr()
    mat=hstack([csr_matrix(G.reshape(-1,ne)),Y,csr_matrix((len(row),1))],format='csr')
    grp=coo_matrix((np.ones(ynum),(np.repeat(np.arange(ng),no),np.arange(ynum))),shape=(ng,ynum)).tocsr()
    mat2=hstack([csr_matrix((ng,ne)),grp,-np.ones((ng,1))],format='csr')
    from scipy.sparse import vstack
    eq=csr_matrix((np.ones(ne),(np.zeros(ne,dtype=int),np.arange(ne))),shape=(1,ne+ynum+1))
    lp=linprog(np.r_[np.zeros(ne+ynum),1.],A_ub=vstack([mat,mat2]),b_ub=np.zeros(mat.shape[0]+ng),A_eq=eq,b_eq=[1.],bounds=[(0,None)]*(ne+ynum+1),method='highs')
    if not lp.success:raise RuntimeError(lp.message)
    w=lp.x[:ne];val=np.max(np.einsum('gobe,e->gob',G,w).max(axis=2).sum(axis=1))
    print('WITNESS',n,val,[(e.tolist(),float(v)) for e,v in zip(errors,w) if v>1e-8],flush=True)
    res={'n':n,'status':'FLOATING_LP_WITNESS','upper_float':float(val),'errors':errors.tolist(),'weights':w.tolist(),'shape':G.shape,'support_count':int(np.sum(w>1e-8))}
    (BASE/'results'/f'character_witness_n{n}.json').write_text(json.dumps(res,indent=2))
    np.savez_compressed(BASE/'results'/f'character_lp_n{n}.npz',G=G,weights=w,labels=np.array([str(x) for x in labels]),errors=errors,dual=lp.ineqlin.marginals)
    return res

if __name__=='__main__':
    for n in map(int,sys.argv[1:] or [4]):search(n)
