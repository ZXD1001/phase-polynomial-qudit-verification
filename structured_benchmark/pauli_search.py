#!/usr/bin/env python3
"""Finite product-Pauli benchmarks. Requires NumPy only; no solver install.

This is a numerical primal/dual search, NOT an interval-certified SDP solver.
The n=2 result has a separate exact arithmetic certificate.
"""
import itertools
import numpy as np

def tensor(matrices):
 result=np.ones((1,1),complex)
 for matrix in matrices:result=np.kron(result,matrix)
 return result

def bases(d):
 if d==2:
  return [np.eye(2),np.array([[1,1],[1,-1]])/np.sqrt(2),np.array([[1,1],[1j,-1j]])/np.sqrt(2)]
 w=np.exp(2j*np.pi/d);x=np.arange(d)[:,None];j=np.arange(d)[None,:]
 return [np.eye(d)]+[w**((t*x*x+j*x)%d)/np.sqrt(d) for t in range(d)]

def tests(d,n):
 labels=np.array(list(itertools.product(range(d),repeat=n)))
 phase=np.sum(labels[:,:-1]**(1 if d==2 else 2)*labels[:,1:],axis=1)%d
 psi=np.exp(2j*np.pi/d*phase)/np.sqrt(d**n)
 Q=np.linalg.qr(np.column_stack([psi,np.eye(d**n)[:,1:]]))[0][:,1:]
 found=[];diagnostics={'maximum_zero_amplitude':0.,'minimum_nonzero_amplitude':1.,'maximum_projector_residual':0.,'maximum_completeness_residual':0.,'maximum_basis_unitarity_residual':0.,'total_bases':(d+1)**n}
 for key in itertools.product(range(d+1),repeat=n):
  basis=tensor([bases(d)[k] for k in key]);amp=basis.conj().T@psi;keep=np.abs(amp)>1e-9
  if np.any(~keep):diagnostics['maximum_zero_amplitude']=max(diagnostics['maximum_zero_amplitude'],float(np.max(np.abs(amp[~keep]))))
  diagnostics['minimum_nonzero_amplitude']=min(diagnostics['minimum_nonzero_amplitude'],float(np.min(np.abs(amp[keep]))))
  P=basis[:,keep]@basis[:,keep].conj().T
  diagnostics['maximum_projector_residual']=max(diagnostics['maximum_projector_residual'],float(np.linalg.norm(P@P-P,2)))
  diagnostics['maximum_completeness_residual']=max(diagnostics['maximum_completeness_residual'],float(np.linalg.norm(P@psi-psi)))
  diagnostics['maximum_basis_unitarity_residual']=max(diagnostics['maximum_basis_unitarity_residual'],float(np.linalg.norm(basis.conj().T@basis-np.eye(d**n),2)))
  if sum(keep)<d**n:found.append((key,P,Q.conj().T@(np.eye(d**n)-P)@Q,int(sum(keep))))
 diagnostics['nontrivial_tests']=len(found)
 assert diagnostics['maximum_zero_amplitude']<1e-12 and diagnostics['minimum_nonzero_amplitude']>1e-4
 assert diagnostics['maximum_projector_residual']<1e-12 and diagnostics['maximum_completeness_residual']<1e-12
 assert np.linalg.norm(Q.conj().T@Q-np.eye(d**n-1))<1e-12
 return psi,Q,found,diagnostics

def simplex(x):
 u=np.sort(x)[::-1];css=np.cumsum(u)-1;ind=np.arange(1,len(u)+1);rho=np.flatnonzero(u-css/ind>0)[-1]
 return np.maximum(x-css[rho]/(rho+1),0)

def smooth(w,R,tau):
 A=np.einsum('i,ijk->jk',w,R,optimize=False);val,vec=np.linalg.eigh(A)
 prob=np.exp(-(val-val[0])/tau);prob/=sum(prob);rho=(vec*prob)@vec.conj().T
 f=float(val[0]-tau*np.log(np.exp(-(val-val[0])/tau).sum()))
 grad=np.einsum('ijk,kj->i',R,rho,optimize=False).real
 return f,grad,float(val[0]),float(max(grad)),rho

def optimize(R,niter=1500):
 # Accelerated projected ascent with backtracking for the concave smoothed
 # minimum eigenvalue. We retain the best primal and dual separately.
 w=np.ones(len(R))/len(R);history=[];bestlower=-1.;bestupper=2.
 for tau in [0.05,0.01,0.002,0.0004,0.00008,0.000016,0.0000032]:
  y=w.copy();theta=1.;L=1.
  for i in range(niter):
   fy,g,*_=smooth(y,R,tau)
   while True:
    wn=simplex(y+g/L);fn,*_=smooth(wn,R,tau);delta=wn-y
    if fn >= fy+g@delta-L*np.dot(delta,delta)/2-1e-13:break
    L*=2
   newtheta=(1+np.sqrt(1+4*theta*theta))/2;y=wn+(theta-1)/newtheta*(wn-w);w=wn;theta=newtheta;L=max(L/1.1,1e-6)
   if i%100==0 or i==niter-1:
    _,_,lower,upper,rho=smooth(w,R,tau)
    if lower>bestlower:bestlower=lower;bestw=w.copy()
    if upper<bestupper:bestupper=upper;bestrho=rho.copy()
  history.append({'tau':tau,'iterations':niter,'best_primal_gap':bestlower,'best_dual_upper':bestupper})
 return bestw,bestrho,history

def complex_json(a):return {'real':a.real.tolist(),'imag':a.imag.tolist()}

def run_n3():
 psi,Q,ss,diag=tests(3,3);R=np.array([s[2] for s in ss]);w,rho,h=optimize(R)
 # Recheck certificates independently of optimizer state.
 w=np.maximum(w,0);w/=sum(w)
 eig,U=np.linalg.eigh((rho+rho.conj().T)/2);eig=np.maximum(eig,0);eig/=sum(eig);rho=(U*eig)@U.conj().T
 Ravg=np.einsum('i,ijk->jk',w,R,optimize=False)
 lower=float(np.linalg.eigvalsh(Ravg)[0]);upper=float(np.max(np.einsum('ijk,kj->i',R,rho,optimize=False).real))
 assert upper-lower<1e-5 and lower>1/6
 # Two-color qubit factor suffices: its gap is 1/2 > qutrit upper bound.
 _,_,qb,qdiag=tests(2,3);qmap={k:r for k,p,r,rank in qb};Rqb=(qmap[(1,0,1)]+qmap[(0,1,0)])/2
 qgap=float(np.linalg.eigvalsh(Rqb)[0]);assert abs(qgap-.5)<1e-12 and upper<qgap
 return {'status':'floating-point primal/dual numerical certificate; not interval-certified',
         'qutrit_diagnostics':diag,'qubit_diagnostics':qdiag,'qubit_two_color_gap':qgap,
         'qutrit_primal_gap':lower,'qutrit_dual_upper':upper,'bracket_width':upper-lower,
         'six_level_full_CRT_class_numeric_bracket':[lower,upper],
         'dual_trace_error':float(abs(np.trace(rho)-1)),'dual_min_eigenvalue':float(np.linalg.eigvalsh(rho)[0]),
         'convergence':h,'weights':[{'basis':list(s[0]),'weight':float(a),'pass_rank':s[3]} for s,a in zip(ss,w)],
         'dual_reduced_density':complex_json(rho),'complement_basis_Q':complex_json(Q),
         'qutrit_target':complex_json(psi),'nonzero_qutrit_setting_count':int(sum(w>1e-12)),
         'construction_six_level_setting_count':int(2*sum(w>1e-12))}
