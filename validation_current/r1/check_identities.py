"""Numerical cross-checks of the analytically proved identities; not proofs.
Reproduce all uniform Z/F comparisons by independent reconstruction, n=2,3,4;
check the n=3 exact-operator claim against a separate native six-level target.
GF(2) rank and compatibility statements are checked with exact binary arithmetic.
"""
if not __debug__:
    raise RuntimeError("Exact checks require normal Python; do not use -O or -OO.")

import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
from itertools import product
import numpy as np
from scipy.linalg import eigh
import json
from lz_optimize import build
OUT=Path(__file__).resolve().parent

def fullgap(Q,a,psi,V,p):
 om=np.einsum('i,ijk->jk',p,Q)
 beta3=eigh(V.conj().T@om@V,eigvals_only=True)[-1]
 beta2=np.max(p@a[:,1:])
 return float(1-max(beta3,beta2))
rows=[]
for n in (2,3,4):
 labels,Q,a,psi,V=build(n);M=len(labels)
 for k in range(n):
  inds=[i for i,(kk,s) in enumerate(labels) if kk==k];p=np.zeros(M);p[inds]=1/len(inds)
  cm=[i for i,(kk,s) in enumerate(labels) if kk==k and len(set(s))==1];pc=np.zeros(M);pc[cm]=1/len(cm)
  rows.append({'n':n,'terminal_0_indexed':k,'MUB':fullgap(Q,a,psi,V,p),'CMUB':fullgap(Q,a,psi,V,pc)})
 p=np.ones(M)/M;cm=[i for i,(kk,s) in enumerate(labels) if len(set(s))==1];pc=np.zeros(M);pc[cm]=1/len(cm)
 rows.append({'n':n,'SMUB':fullgap(Q,a,psi,V,p),'SCMUB':fullgap(Q,a,psi,V,pc)})
# Independent native d=6 construction at n=3, no CRT used below.
d=6;n=3;w=np.exp(2j*np.pi/6);X=np.array(list(product(range(d),repeat=n)))
g=lambda x:(2*x*x+3*x)%6
phase=(g(X[:,0])*X[:,1]+g(X[:,1])*X[:,2])%6
psi=w**phase/d**(n/2)
I=np.eye(d);F=w**np.outer(np.arange(d),np.arange(d))/np.sqrt(d)
def kr(mats):
 z=np.array([[1.]])
 for m in mats:z=np.kron(z,m)
 return z
Pvs=[]
for k in range(n):
 mats=[I]*n;mats[k]=np.ones((d,d))/d;pp=kr(mats)
 Pvs.append(w**phase[:,None]*pp*w**(-phase[None,:]))
def nativeQ(k,strings):
 B=[i for i in range(n) if i!=k];mats=[I]*n
 for b,s in zip(B,strings):mats[b]=I if s==0 else F
 U=kr(mats);trans=(U.conj().T@psi).reshape((d,)*n);C=[]
 for o in product(range(d),repeat=n-1):
  sl=[slice(None) if i==k else o[B.index(i)] for i in range(n)]
  v=trans[tuple(sl)]
  if np.linalg.norm(v)<1e-12:continue
  v=v/np.linalg.norm(v);col=np.zeros((d,)*n,complex);col[tuple(sl)]=v;C.append(U@col.ravel())
 C=np.array(C).T;return C@C.conj().T
q0=nativeQ(0,(0,1));q1=nativeQ(1,(0,0))
res0=float(np.linalg.norm(q0-Pvs[0]@Pvs[2],2));res1=float(np.linalg.norm(q1-Pvs[1],2));
assert res0<1e-12 and res1<1e-12

def rank2(A):
 A=np.array(A,dtype=int).copy()%2;r=0
 for c in range(A.shape[1]):
  nz=np.flatnonzero(A[r:,c]);
  if len(nz)==0:continue
  p=r+nz[0];A[[r,p]]=A[[p,r]]
  for j in range(A.shape[0]):
   if j!=r and A[j,c]:A[j]^=A[r]
  r+=1
  if r==A.shape[0]:break
 return r
ranks=[]
for n in range(2,31):
 A=np.zeros((n,n),int)
 for j in range(n-1):A[j,j+1]=A[j+1,j]=1
 for k in range(n):
  r=rank2(np.delete(A,k,axis=0));expect=n-1 if n%2==0 or k%2==0 else n-2
  assert r==expect
  ranks.append({'n':n,'k':k,'rank':r,'common_qubit_dimension':2**(r-1)})
# finite checks of the all-n 1/2 upper bound: compatibility subspace does not
# reject both Z0|G> and Z1|G> within any one canonical real-front test.
for n in range(2,8):
 for k in range(n):
  B=[j for j in range(n) if j!=k]
  for s in product(range(2),repeat=n-1):
   detected=[False,False]
   for aa in product(range(2),repeat=n):
    bb=[((aa[j-1] if j else 0)+(aa[j+1] if j<n-1 else 0))%2 for j in range(n)]
    if any(aa[j] if t==0 else bb[j] for j,t in zip(B,s)):continue
    assert all(aa[j]*bb[j]==0 for j in range(n))
    detected[0]|=bool(aa[0]);detected[1]|=bool(aa[1])
   assert not all(detected)
res={'status':'PASS','uniform_rows':rows,'native_n3_residual_Q0ZF_minus_P0P2':res0,'native_n3_residual_Q1ZZ_minus_P1':res1,'rank_formula_exact_finite_checks_n2_to_30':ranks,'half_upper_compatibility_finite_checks_n2_to_7':True}
(OUT/'identity_checks.json').write_text(json.dumps(res,indent=2))
print(json.dumps({k:v for k,v in res.items() if k!='rank_formula_exact_finite_checks_n2_to_30'},indent=2))
