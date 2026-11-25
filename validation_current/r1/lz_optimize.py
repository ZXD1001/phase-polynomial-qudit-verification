"""Independent Z/F Li--Zhu reconstruction and weight optimization.
Qutrit branch support: exact integer cyclotomic sums. Optimization: floating LP
outer approximation using minimum-eigenvector cuts (NOT interval certification).
All tested front strings and terminals can be mixed with arbitrary weights.
A commuting qubit stabilizer factor reduces the eigenproblems to dimension 3^n.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from itertools import product
from pathlib import Path
import numpy as np
from scipy.linalg import eigh, null_space
from scipy.optimize import linprog
import json,sys,time
OUT=Path(__file__).resolve().parent

def configs(d,n):return np.array(list(product(range(d),repeat=n)),dtype=int)
def target_phase(X,d):
    if d==2:return np.sum(X[:,:-1]*X[:,1:],axis=1)%2
    return np.sum(X[:,:-1]**2*X[:,1:],axis=1)%3

def conditional_test(d,n,k,s):
    X=configs(d,n); f=target_phase(X,d); root=np.exp(2j*np.pi/d)
    B=[i for i in range(n) if i!=k]; s=np.array(s); nF=sum(s)
    cols=[];support=[]
    for o in product(range(d),repeat=n-1):
        o=np.array(o);allow=np.ones(len(X),dtype=bool);phase=np.zeros(len(X),dtype=int)
        for j,site in enumerate(B):
            if s[j]==0:allow&=(X[:,site]==o[j])
            else:phase+=o[j]*X[:,site]
        exps=(f-phase)%d
        coeff=[]
        for y in range(d):
            counts=np.bincount(exps[allow & (X[:,k]==y)],minlength=d)
            if d==2: coeff.append((int(counts[0]-counts[1]),0))
            else: coeff.append((int(counts[0]-counts[2]),int(counts[1]-counts[2])))
        if all(a==b==0 for a,b in coeff):continue
        v=np.array([a+b*root for a,b in coeff]);v/=np.linalg.norm(v)
        col=np.zeros(d**n,dtype=complex)
        col[allow]=root**(phase[allow]%d)*v[X[allow,k]]/d**(nF/2)
        cols.append(col);support.append({'outcome':o.tolist(),'coeff':coeff})
    C=np.array(cols).T
    Q=C@C.conj().T
    psi=root**f/d**(n/2)
    assert np.linalg.norm(Q@psi-psi)<1e-10
    assert np.linalg.norm(C.conj().T@C-np.eye(C.shape[1]))<1e-10
    return Q,support

def build(n):
    cache=OUT/f'lz_inputs_n{n}.npz'
    labels=[(k,s) for k in range(n) for s in product(range(2),repeat=n-1)]
    if cache.exists():
        dat=np.load(cache);return labels,dat['Q3'],dat['a2'],dat['psi3'],dat['V']
    X2=configs(2,n);f2=target_phase(X2,2)
    G=(-1.)**f2[:,None]*(-1.)**(X2@X2.T)/2**(n/2)
    Q3=[];a2=[];branches=[];resid=0
    for k,s in labels:
        q2,_=conditional_test(2,n,k,s);q3,sp=conditional_test(3,n,k,s)
        q2sector=G.T@q2@G;vals=np.diag(q2sector).real
        resid=max(resid,float(np.linalg.norm(q2sector-np.diag(np.rint(vals)))))
        assert resid<1e-10
        a2.append(np.rint(vals).astype(int));Q3.append(q3);branches.append(sp)
    X3=configs(3,n);psi3=np.exp(2j*np.pi*target_phase(X3,3)/3)/3**(n/2)
    V=null_space(psi3.conj()[None,:])
    np.savez_compressed(cache,Q3=np.array(Q3),a2=np.array(a2),psi3=psi3,V=V)
    (OUT/f'lz_branches_n{n}.json').write_text(json.dumps({'labels':labels,'qutrit_branches':branches,'qubit_diagonalization_residual':resid}))
    return labels,np.array(Q3),np.array(a2),psi3,V

def optimize(n,indices=None,maxit=1200,tol=2e-7):
    labels,Q3,a2,psi,V=build(n)
    if indices is None:indices=list(range(len(labels)))
    Q3=Q3[indices];a2=a2[indices];labels=[labels[i] for i in indices]
    M=len(labels);D=len(psi)-1
    r=np.array([V.conj().T@(np.eye(len(psi))-q)@V for q in Q3]); r=(r+r.conj().transpose(0,2,1))/2
    rej2=1-a2[:,1:].T
    cuts=list(rej2.astype(float));cut_vectors=[None]*len(cuts)
    p=np.ones(M)/M;bestp=p.copy();low=-1;history=[]
    for it in range(maxit):
        R=np.einsum('i,ijk->jk',p,r)
        vals,vecs=eigh(R,subset_by_index=(0,min(7,D-1)),driver='evr')
        val=min(float(vals[0]),float(np.min(rej2@p)))
        if val>low:low=val;bestp=p.copy()
        # Add all almost-low directions; full few eigenvectors helps degeneracy.
        for j in range(vecs.shape[1]):
            u=vecs[:,j]
            z=np.einsum('j,ijk,k->i',u.conj(),r,u).real
            cuts.append(z);cut_vectors.append(u.copy())
        C=np.array(cuts)
        lp=linprog(np.r_[np.zeros(M),-1.],A_ub=np.c_[-C,np.ones(len(C))],b_ub=np.zeros(len(C)),A_eq=np.array([np.r_[np.ones(M),0.]]),b_eq=[1.],bounds=[(0,1)]*M+[(0,1)],method='highs')
        if not lp.success:raise RuntimeError(lp.message)
        p=np.maximum(lp.x[:-1],0);p/=p.sum();upper=float(lp.x[-1])
        history.append([it,low,upper])
        if it%25==0:print('n',n,'iter',it,'lower',low,'outer_upper',upper,flush=True)
        if upper-low<tol:break
    # Export actual PSD dual from nonnegative LP multipliers, split across blocks.
    lamb=np.maximum(-lp.ineqlin.marginals,0);lamb/=lamb.sum()
    dual3=np.zeros((D,D),complex);dual2=np.zeros(rej2.shape[0]);
    for j,(w,u) in enumerate(zip(lamb,cut_vectors)):
        if u is None:dual2[j]=w
        else:dual3+=w*np.outer(u,u.conj())
    exp=np.einsum('ijk,kj->i',r,dual3).real+rej2.T@dual2
    dual_upper=float(exp.max())
    result={'n':n,'status':'FLOATING_PRIMAL_DUAL_NOT_INTERVAL_CERTIFIED','labels':labels,'weights':bestp.tolist(),'primal_gap':low,'dual_upper':dual_upper,'iterations':len(history),'dual3_trace':float(np.trace(dual3).real),'dual2_weight':float(dual2.sum()),'dual3_min_eigenvalue':float(eigh(dual3,eigvals_only=True)[0]),'history':history}
    name='all' if len(indices)==n*2**(n-1) else '_'.join(map(str,indices))
    (OUT/f'lz_optimized_n{n}_{name}.json').write_text(json.dumps(result,indent=2))
    np.savez_compressed(OUT/f'lz_dual_n{n}_{name}.npz',dual3=dual3,dual2=dual2,V=V,expectations=exp)
    print('FINAL',n,low,dual_upper,'support',[(labels[i],round(v,9)) for i,v in enumerate(bestp) if v>1e-6],flush=True)
    return result

if __name__=='__main__':
    for n in map(int,sys.argv[1:] or [2,3,4]):optimize(n)
