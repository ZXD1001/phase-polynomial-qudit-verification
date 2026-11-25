#!/usr/bin/env python3
"""Full bounded reproduction entry. Python + NumPy only; writes beside itself.
Original manuscripts, prior reports and Zotero are never read or modified.
"""
import sys
sys.dont_write_bytecode=True
import itertools,json,platform,time
from pathlib import Path
import numpy as np
from exact_n2_certificate import certificate
from pauli_search import tensor,tests,run_n3

W=np.exp(2j*np.pi/6);PLUS=np.ones(6,complex)/np.sqrt(6);PPLUS=np.outer(PLUS,PLUS)
def g(x):return (2*x*x+3*x)%6

def physical(n,A,setting=None):
 """Construct actual measurement effect from conditional outcome rays.
 setting=None: adaptive. Dict v:t: preset local menus, matching-only accept.
 A is independent; complement Z first, then local rays (or ignore mismatch).
 Terminal n-1 has its actual one-element menu t=0.
 """
 B=[v for v in range(n) if v not in A];D=6**n;out=np.zeros((D,D),complex)
 labels=np.array(list(itertools.product(range(6),repeat=n)))
 for record in itertools.product(range(6),repeat=len(B)):
  m=dict(zip(B,record));idx=np.flatnonzero(np.all(labels[:,B]==np.array(record),axis=1))
  factors=[]
  for v in A:
   left=m.get(v-1,0);right=m.get(v+1,0);required=right%3
   t=required if setting is None else setting[v]
   if t!=required:factors.append(np.eye(6));continue
   j=(-(g(left)+3*right))%6;x=np.arange(6)
   ray=W**((2*t*x*x-j*x)%6)/np.sqrt(6);factors.append(np.outer(ray,ray.conj()))
  out[np.ix_(idx,idx)]=tensor(factors)
 return out

def path_check(n):
 labels=np.array(list(itertools.product(range(6),repeat=n)))
 phase=np.sum(g(labels[:,:-1])*labels[:,1:],axis=1)%6;u=W**phase;D=6**n
 def dressed(M):return u[:,None]*M*u.conj()[None,:]
 ad=[];ind=[];shared=[];err_ad=[];err_ind=[];settings=[]
 for parity in (0,1):
  A=list(range(parity,n,2));q={v:1. if v==n-1 else 1/3 for v in A}
  ideal=dressed(tensor([PPLUS if v in A else np.eye(6) for v in range(n)]));actual=physical(n,A)
  err_ad.append(float(np.max(np.abs(actual-ideal))));ad.append(ideal)
  menus=[(0,) if v==n-1 else (0,1,2) for v in A];combos=list(itertools.product(*menus));settings.append(len(combos))
  physical_average=np.zeros((D,D),complex)
  for combo in combos:physical_average+=physical(n,A,dict(zip(A,combo)))/len(combos)
  formula=dressed(tensor([(1-q[v])*np.eye(6)+q[v]*PPLUS if v in A else np.eye(6) for v in range(n)]))
  err_ind.append(float(np.max(np.abs(physical_average-formula))));ind.append(formula)
  shared.append(sum(physical(n,A,{v:0 if v==n-1 else t for v in A}) for t in range(3))/3)
 omega_ad=sum(ad)/2;omega_ind=sum(ind)/2;omega_shared=sum(shared)/2
 eta=u*W**(2*labels[:,0])/np.sqrt(D)
 sat_ind=float(np.max(np.abs(omega_ind@eta-(5/6)*eta)));sat_shared=float(np.max(np.abs(omega_shared@eta-(5/6)*eta)))
 result={'n':n,'hilbert_dimension':D,'actual_independent_setting_counts_by_color':settings,
         'independent_sampling_random_trits_per_copy_at_most':max(sum(v<n-1 for v in range(parity,n,2)) for parity in (0,1)),
         'adaptive_physical_vs_operator_max_errors':err_ad,'independent_physical_vs_HT_product_formula_max_errors':err_ind,
         'independent_vs_shared_max_entry_difference':float(np.max(np.abs(omega_ind-omega_shared))),
         'independent_saturating_vector_residual':sat_ind,'shared_saturating_vector_residual':sat_shared,
         'analytical_adaptive_gap':'1/2','analytical_uniform_color_independent_gap':'1/6','analytical_uniform_color_shared_gap':'1/6'}
 assert max(err_ad+err_ind+[sat_ind,sat_shared])<1e-12
 if n<=3:
  def gap(M):return float(1-np.linalg.eigvalsh(M)[-2])
  result.update(adaptive_gap_eigensolver=gap(omega_ad),independent_gap_eigensolver=gap(omega_ind),shared_gap_eigensolver=gap(omega_shared))
  assert abs(result['adaptive_gap_eigensolver']-.5)<1e-12 and abs(result['independent_gap_eigensolver']-1/6)<1e-12 and abs(result['shared_gap_eigensolver']-1/6)<1e-12
  if n==2:
   result['optimized_four_setting_gap']=gap(.75*ind[0]+.25*ind[1]);assert abs(result['optimized_four_setting_gap']-.25)<1e-12
 return result

def full_crt_n2_check():
 # Direct physical-six-level construction checks all 144 product bases.
 qb,Qb,sb,db=tests(2,2);qt,Qt,st,dt=tests(3,2)
 from pauli_search import bases
 target=np.array([W**(g(x)*y%6) for x,y in itertools.product(range(6),repeat=2)])/6
 crt=np.zeros((36,36))
 # CRT order is (a0,a1,b0,b1); physical order is (x0,x1).
 for x,y in itertools.product(range(6),repeat=2):crt[6*x+y,(2*(x%2)+y%2)*9+3*(x%3)+y%3]=1
 assert np.max(np.abs(target-crt@np.kron(qb,qt)))<1e-12
 max_tensor_error=0.;omega=np.zeros((36,36),complex);weights_sum=0.;count=0
 def canonical(psi,B):
  keep=np.abs(B.conj().T@psi)>1e-9;return B[:,keep]@B[:,keep].conj().T
 for kb in itertools.product(range(3),repeat=2):
  Bb=tensor([bases(2)[i] for i in kb]);Pb=canonical(qb,Bb)
  wb=1/3 if kb in [(0,1),(1,0),(2,2)] else 0.
  for kt in itertools.product(range(4),repeat=2):
   Bt=tensor([bases(3)[i] for i in kt]);Pt=canonical(qt,Bt)
   Bphysical=crt@np.kron(Bb,Bt);Pphysical=canonical(target,Bphysical);Pfactor=crt@np.kron(Pb,Pt)@crt.T
   max_tensor_error=max(max_tensor_error,float(np.max(np.abs(Pphysical-Pfactor))))
   wt=2/11 if kt==(0,1) else 1/11 if kt[0] in (1,2,3) and kt[1] in (0,2,3) else 0.
   if wb*wt:omega+=wb*wt*Pphysical;weights_sum+=wb*wt;count+=1
 gap=float(1-np.linalg.eigvalsh(omega)[-2]);assert abs(gap-3/11)<1e-12 and max_tensor_error<1e-12
 return {'enumerated_physical_six_level_bases':144,'max_direct_physical_vs_CRT_tensor_acceptance_error':max_tensor_error,
         'construction_settings':count,'weight_sum':weights_sum,'direct_physical_gap':gap,
         'exact_optimum_gap_in_full_defined_CRT_class':'3/11','qubit_diagnostics':db,'qutrit_diagnostics':dt}

if __name__=='__main__':
 start=time.perf_counter();out=Path(__file__).resolve().parent
 result={'python':platform.python_version(),'numpy':np.__version__,'scope':'six-level path; matching and fixed local CRT product-Pauli finite classes',
         'paths':[path_check(n) for n in (2,3,4)],'n2_exact_qutrit_certificate':certificate(),
         'n2_direct_six_level':full_crt_n2_check(),'n3_numeric':run_n3()}
 result['runtime_seconds']=time.perf_counter()-start
 (out/'benchmark_results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
 summary={k:result[k] for k in ['python','numpy','paths','n2_direct_six_level','runtime_seconds']}
 summary['n3_numeric_bracket']=result['n3_numeric']['six_level_full_CRT_class_numeric_bracket']
 summary['n3_numerical_bracket_width']=result['n3_numeric']['bracket_width']
 print(json.dumps(summary,ensure_ascii=False,indent=2))
