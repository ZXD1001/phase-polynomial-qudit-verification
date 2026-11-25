#!/usr/bin/env python3
"""Exact Q(omega_3) arithmetic certificate; Python standard library only.

All equalities, traces, and projector conditions here are rational identities,
not tolerance-based floating point tests. Basis 0 is Z, bases 1+t have columns
omega^(t*x*x+j*x)/sqrt(3). Lexicographic computational order x then y.
"""
from fractions import Fraction as F
from itertools import product
import json

Z=(F(0),F(0)); ONE=(F(1),F(0)); W=(F(0),F(1))
def add(a,b):return (a[0]+b[0],a[1]+b[1])
def neg(a):return (-a[0],-a[1])
def mul(a,b):return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]-a[1]*b[1])
def scale(a,q):return (a[0]*q,a[1]*q)
def wp(k):return [ONE,W,(F(-1),F(-1))][k%3]
def zero(n=9):return [[Z for _ in range(n)] for _ in range(n)]
def eye(n=9):return [[ONE if i==j else Z for j in range(n)] for i in range(n)]
def plus(A,B):return [[add(a,b) for a,b in zip(x,y)] for x,y in zip(A,B)]
def scal(A,q):return [[scale(a,q) for a in x] for x in A]
def sub(A,B):return plus(A,scal(B,F(-1)))
def mm(A,B):
 C=zero(len(A))
 for i in range(len(A)):
  for j in range(len(A)):
   for k in range(len(A)):C[i][j]=add(C[i][j],mul(A[i][k],B[k][j]))
 return C
def tr(A):
 a=Z
 for i in range(len(A)):a=add(a,A[i][i])
 assert a[1]==0
 return a[0]
def power(A,k):
 B=eye(len(A))
 for _ in range(k):B=mm(B,A)
 return B

def certificate():
 labels=list(product(range(3),repeat=2));phase=[x*x*y%3 for x,y in labels]
 psi=[[scale(wp(phase[i]-phase[j]),F(1,9)) for j in range(9)] for i in range(9)]
 tests={};ranks={};I=eye()
 for key in product(range(4),repeat=2):
  P=zero();rank=0;k=sum(t!=0 for t in key)
  for js in labels:
   vec=[]
   for xs in labels:
    allowed=all(t!=0 or x==j for t,x,j in zip(key,xs,js))
    phase_b=sum((t-1)*x*x+j*x for t,x,j in zip(key,xs,js) if t!=0)%3
    vec.append(phase_b if allowed else None)
   amplitude=Z
   for i,b in enumerate(vec):
    if b is not None:amplitude=add(amplitude,wp(phase[i]-b))
   if amplitude==Z:continue
   rank+=1
   for i,a in enumerate(vec):
    for j,b in enumerate(vec):
     if a is not None and b is not None:P[i][j]=add(P[i][j],scale(wp(a-b),F(1,3**k)))
  assert mm(P,P)==P and mm(P,psi)==psi
  tests[key]=sub(I,P);ranks[key]=rank
 weights={k:(F(2,11) if k==(0,1) else F(1,11) if k[0] in (1,2,3) and k[1] in (0,2,3) else F(0)) for k in tests}
 R=zero()
 for k in tests:R=plus(R,scal(tests[k],weights[k]))
 assert sum(weights.values())==1
 # Exact characteristic annihilator, followed by unique-zero eigenspace.
 factors=[sub(R,scal(I,F(k,11))) for k in (0,3,4,5)]
 poly=eye()
 for a in factors:poly=mm(poly,a)
 assert poly==zero()
 E0=eye()
 for k in (3,4,5):E0=mm(E0,scal(sub(R,scal(I,F(k,11))),F(-11,k)))
 assert E0==psi
 # Minimum-eigenvalue projector E, and its commuting rank-two part EL.
 E=sub(scal(mm(factors[2],factors[3]),F(121,2)),scal(psi,F(10)))
 L=tests[(0,1)];EL=mm(E,L)
 assert mm(E,E)==E and tr(E)==4
 assert mm(L,E)==EL and mm(EL,EL)==EL and tr(EL)==2
 sigma=sub(scal(E,F(4,11)),scal(EL,F(5,22)))
 assert tr(sigma)==1 and mm(psi,sigma)==zero()
 # sigma=4/11(E-EL)+3/22 EL is positive and supported off the target.
 vals={k:tr(mm(tests[k],sigma)) for k in tests}
 assert max(vals.values())==F(3,11)
 for k,v in vals.items():
  expected=F(3,11) if weights[k]>0 else F(3,22) if k[0] in (1,2,3) and k[1]==1 else F(0)
  assert v==expected
 return {'arithmetic':'exact fractions in Q(omega_3)', 'all_16_tests_projective_and_perfectly_complete':True,
         'primal_weights':{str(k):str(v) for k,v in weights.items() if v},
         'annihilating_spectrum':['0','3/11','4/11','5/11'], 'kernel_projector_equals_target':True,
         'minimum_spectral_projector_rank':int(tr(E)), 'dual_rank_two_projector_rank':int(tr(EL)),
         'dual_eigenvalues_nonzero':['4/11','4/11','3/22','3/22'],
         'all_dual_rejection_expectations':{str(k):str(v) for k,v in vals.items()},
         'certified_optimum_gap':'3/11','certified_class':'all 4^2 fixed product qutrit Pauli bases, arbitrary classical perfect-completeness acceptance and setting probabilities'}

if __name__=='__main__':print(json.dumps(certificate(),indent=2))
