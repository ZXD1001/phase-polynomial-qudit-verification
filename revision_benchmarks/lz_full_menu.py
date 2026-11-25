"""Exploratory canonical conditional tests with front menu Z,B0,B1,B2.
Each native-six-level front basis is the manuscript's calibrated local PVM.
CRT: Z2*Z3 or X2*quadratic-chirp qutrit Pauli basis. Floating optimization.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np
from scipy.linalg import null_space,eigh
from scipy.optimize import linprog
from itertools import product
from pathlib import Path
import json,sys,time
from lz_optimize import configs,target_phase,conditional_test
OUT=Path(__file__).resolve().parent

def q3test(n,k,s):
 X=configs(3,n);f=target_phase(X,3);root=np.exp(2j*np.pi/3)
 B=[i for i in range(n) if i!=k];nf=sum(t!=0 for t in s);cols=[]
 for o in product(range(3),repeat=n-1):
  allow=np.ones(len(X),bool);phase=np.zeros(len(X),int)
  for t,site,oj in zip(s,B,o):
   if not t:allow&=X[:,site]==oj
   else:phase+=(t-1)*X[:,site]**2+oj*X[:,site]
  e=(f-phase)%3;coeff=[]
  for y in range(3):
   c=np.bincount(e[allow&(X[:,k]==y)],minlength=3)
   coeff.append((int(c[0]-c[2]),int(c[1]-c[2])))
  if all(a==b==0 for a,b in coeff):continue
  v=np.array([a+b*root for a,b in coeff]);v/=np.linalg.norm(v)
  col=np.zeros(len(X),complex);col[allow]=root**(phase[allow]%3)*v[X[allow,k]]/3**(nf/2);cols.append(col)
 C=np.array(cols).T;Q=C@C.conj().T;psi=root**f/3**(n/2)
 assert np.linalg.norm(Q@psi-psi)<1e-10
 return Q

def build(n):
 labels=[(k,s) for k in range(n) for s in product(range(4),repeat=n-1)]
 cache=OUT/f'lz_fullmenu_inputs_n{n}.npz'
 if cache.exists():
  data=np.load(cache);return labels,data['Q3'],data['a2'],data['V'],data['psi3']
 X2=configs(2,n);f2=target_phase(X2,2);G=(-1.)**f2[:,None]*(-1.)**(X2@X2.T)/2**(n/2)
 a2=[];Q3=[];c2={}
 for k,s in labels:
  s2=tuple(int(t!=0) for t in s)
  if (k,s2) not in c2:
   Q,_=conditional_test(2,n,k,s2);M=G.T@Q@G;v=np.rint(np.diag(M).real).astype(int)
   assert np.linalg.norm(M-np.diag(v))<1e-10;c2[k,s2]=v
  a2.append(c2[k,s2]);Q3.append(q3test(n,k,s))
 X3=configs(3,n);psi3=np.exp(2j*np.pi*target_phase(X3,3)/3)/3**(n/2);V=null_space(psi3.conj()[None,:])
 Q3=np.array(Q3);a2=np.array(a2);np.savez_compressed(cache,Q3=Q3,a2=a2,V=V,psi3=psi3)
 return labels,Q3,a2,V,psi3

def run(n):
 labels,Q3,a2,V,psi=build(n);D=len(psi)-1;M=len(labels)
 rr=np.array([V.conj().T@(np.eye(len(psi))-q)@V for q in Q3]);rr=(rr+rr.conj().transpose(0,2,1))/2
 rej2=1-a2[:,1:].T;cuts=list(rej2.astype(float));vectors=[None]*len(cuts)
 p=np.ones(M)/M;bestp=p.copy();lo=-1;history=[]
 for it in range(161):
  R=np.einsum('i,ijk->jk',p,rr);vals,vecs=eigh(R,subset_by_index=(0,min(7,D-1)))
  val=min(float(vals[0]),float(np.min(rej2@p)))
  if val>lo:lo=val;bestp=p.copy()
  for u in vecs.T:
   cuts.append(np.einsum('j,ijk,k->i',u.conj(),rr,u).real);vectors.append(u.copy())
  C=np.array(cuts);lp=linprog(np.r_[np.zeros(M),-1.],A_ub=np.c_[-C,np.ones(len(C))],b_ub=np.zeros(len(C)),A_eq=np.array([np.r_[np.ones(M),0.]]),b_eq=[1.],bounds=[(0,1)]*(M+1),method='highs')
  if not lp.success:raise RuntimeError(lp.message)
  p=np.maximum(lp.x[:-1],0);p/=p.sum();up=float(lp.x[-1]);history.append([it,lo,up])
  if it%10==0:
   print(n,it,lo,up,flush=True)
   (OUT/f'lz_fullmenu_checkpoint_n{n}.json').write_text(json.dumps({'n':n,'status':'FLOATING_FEASIBLE_CHECKPOINT_NOT_CERTIFIED','iteration':it,'primal_gap':lo,'outer_upper':up,'labels':labels,'weights':bestp.tolist(),'history':history},indent=2))
  if len(cuts)>450 and up-lo>=5e-7 and it<160:
   active=np.flatnonzero(-lp.ineqlin.marginals>1e-9).tolist()
   keep=sorted(set(list(range(len(rej2)))+active+list(range(max(len(rej2),len(cuts)-120),len(cuts)))))
   cuts=[cuts[j] for j in keep];vectors=[vectors[j] for j in keep]
  if up-lo<5e-7:break
 lam=np.maximum(-lp.ineqlin.marginals,0);lam/=sum(lam);dual3=np.zeros((D,D),complex);dual2=np.zeros(rej2.shape[0])
 for j,(w,u) in enumerate(zip(lam,vectors)):
  if u is None:dual2[j]=w
  else:dual3+=w*np.outer(u,u.conj())
 dualup=float((np.einsum('ijk,kj->i',rr,dual3).real+rej2.T@dual2).max())
 result={'n':n,'status':'FLOATING_PRIMAL_DUAL_NOT_CERTIFIED','menu':'Z,B0,B1,B2; all front assignments; arbitrary terminal/setting weights','primal_gap':lo,'dual_upper':dualup,'labels':labels,'weights':bestp.tolist(),'history':history}
 (OUT/f'lz_fullmenu_n{n}.json').write_text(json.dumps(result,indent=2));np.savez_compressed(OUT/f'lz_fullmenu_dual_n{n}.npz',dual3=dual3,dual2=dual2,V=V)
 print('FINAL',n,lo,dualup,[(labels[i],float(w)) for i,w in enumerate(bestp) if w>1e-7],flush=True)
if __name__=='__main__':
 for n in map(int,sys.argv[1:] or [4]):run(n)
