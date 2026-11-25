"""Produce and independently check an EXACT rational dual upper bound for n=4.
The 32 allowed tests are all terminals and all fixed Z/F front strings.
The candidate is a convex combination of exactly target-orthogonal rays over
Q(w3) in qubit target sector, plus qubit error sectors tensored with H4.
Positivity follows by construction. All expectations use integer/Fraction
arithmetic. Float eigenvectors are used ONLY to propose the input certificate.
With --verify-only, no eigensolver is used.
"""
if not __debug__:
    raise RuntimeError("Exact checks require normal Python; do not use -O or -OO.")

from pathlib import Path
from itertools import product
from fractions import Fraction as F
import json,sys
P=Path(__file__).resolve().parent
ROOT=((1,0),(0,1),(-1,-1))
def add(x,y):return x[0]+y[0],x[1]+y[1]
def neg(x):return -x[0],-x[1]
def mul(x,y):
 a,b=x;c,d=y
 return a*c-b*d,a*d+b*c-b*d
def conj(x):return x[0]-x[1],-x[1]
def norm(x):a,b=x;return a*a-a*b+b*b
def inner(v,w):
 a=b=0
 for x,y in zip(v,w):
  c,d=mul(conj(x),y);a+=c;b+=d
 return a,b
X=list(product(range(3),repeat=4));D=len(X)
h=[ROOT[sum(x[i]**2*x[i+1] for i in range(3))%3] for x in X]
cert_path=P/'lz_n4_exact_upper_certificate.json'
if '--verify-only' not in sys.argv:
 import numpy as np
 from scipy.linalg import eigh
 dat=np.load(P/'lz_dual_n4_all.npz');V=dat['V'];rho=V@dat['dual3']@V.conj().T
 vals,vecs=eigh((rho+rho.conj().T)/2)
 vecsout=[];weights=[];scale=10**7
 for val,v in zip(vals,vecs.T):
  if val<1e-11:continue
  bb=np.rint(v.imag*2/np.sqrt(3)*scale).astype(np.int64)
  aa=np.rint((v.real+v.imag/np.sqrt(3))*scale).astype(np.int64)
  z=[(int(a),int(b)) for a,b in zip(aa,bb)]
  c=inner(h,z)
  q=[add((D*a,D*b),neg(mul(hj,c))) for (a,b),hj in zip(z,h)]
  assert inner(h,q)==(0,0)
  vecsout.append(q);weights.append(F(float(val)).limit_denominator(10**10))
 d2=[]
 for j,val in enumerate(dat['dual2']):
  if val>1e-11:d2.append((j+1,F(float(val)).limit_denominator(10**10)))
 total=sum(weights)+sum(w for j,w in d2)
 c={'n':4,'candidate_generation':'Floating eigenvectors rounded to Q(w3), followed by EXACT target projection. Generation is not certification.',
 'qubit_target_sector_qutrit_vectors':[{'weight':str(w/total),'v':v} for w,v in zip(weights,vecsout)],
 'qubit_error_sectors':[{'sector':j,'weight':str(w/total)} for j,w in d2],
 'claimed_upper':'197537/500000','class':'All 32 canonical minimal Li-Zhu conditional-ray tests, fixed Z/F front bases, arbitrary joint weights over terminal and basis string'}
 cert_path.write_text(json.dumps(c,indent=2))
c=json.loads(cert_path.read_text());vecs=[[(int(a),int(b)) for a,b in row['v']] for row in c['qubit_target_sector_qutrit_vectors']]
w3=[F(row['weight']) for row in c['qubit_target_sector_qutrit_vectors']];d2=[(row['sector'],F(row['weight'])) for row in c['qubit_error_sectors']]
assert sum(w3)+sum(w for _,w in d2)==1 and all(w>=0 for w in w3) and all(w>=0 for _,w in d2)
Nv=[sum(norm(z) for z in v) for v in vecs]
assert all(N>0 for N in Nv) and all(inner(h,v)==(0,0) for v in vecs)

def qubit_accepts(k,s,sector):
 """Exact graph-stabilizer compatibility; no numerical diagonalization."""
 B=[j for j in range(4) if j!=k];err=tuple((sector>>(3-j))&1 for j in range(4))
 for a in product(range(2),repeat=4):
  b=tuple(((a[j-1] if j else 0)+(a[j+1] if j<3 else 0))%2 for j in range(4))
  if any((a[j] if t==0 else b[j]) for j,t in zip(B,s)):continue
  if sum(x*y for x,y in zip(a,err))%2:return 0
 return 1

def exact_rejection(k,s,v,Nv):
 B=[j for j in range(4) if j!=k];nF=sum(s);acc=F(0)
 for o in product(range(3),repeat=3):
  phases=[];q=[(0,0)]*3
  for ix,x in enumerate(X):
   if any(t==0 and x[j]!=oj for j,t,oj in zip(B,s,o)):
    phases.append(None);continue
   ph=sum(oj*x[j] for j,t,oj in zip(B,s,o) if t)%3
   phases.append(ph)
   q[x[k]]=add(q[x[k]],mul(h[ix],ROOT[-ph%3]))
  Nq=sum(norm(z) for z in q)
  if Nq==0:continue
  # |b_front> tensor q / sqrt(3^nF * Nq)
  raw=[(0,0) if ph is None else mul(ROOT[ph],q[x[k]]) for ph,x in zip(phases,X)]
  z=inner(raw,v)
  acc+=F(norm(z),3**nF*Nq*Nv)
 assert 0<=acc<=1
 return 1-acc

expect=[]
for k in range(4):
 for ss in product(range(2),repeat=3):
  e=sum(w*exact_rejection(k,ss,v,nv) for w,v,nv in zip(w3,vecs,Nv))
  e+=sum(w*(1-qubit_accepts(k,ss,j)) for j,w in d2)
  assert e<=F(c['claimed_upper']),(k,ss,float(e))
  expect.append({'terminal':k,'front_string':ss,'expectation_exact':str(e),'decimal_display':float(e)})
mx=max(F(e['expectation_exact']) for e in expect)
res={'status':'EXACT_PASS','upper_bound':'197537/500000','upper_bound_decimal':0.395074,'maximum_exact_expectation':str(mx),'maximum_decimal_display':float(mx),'all_32_expectations':expect,'qubit_target_qutrit_rays':len(vecs),'qubit_error_terms':len(d2),'PSD_reason':'Convex mixture of normalized rays and sector states; all coefficients nonnegative rational. Exact projection enforces target orthogonality.', 'verification_arithmetic':'Python arbitrary-precision integers and Fraction only; eigensolver not used by verification.'}
(P/'lz_n4_exact_upper_result.json').write_text(json.dumps(res,indent=2))
print('EXACT_PASS: all32 expectations <=197537/500000; max display',float(mx),'raycount',len(vecs),'errorsectorcount',len(d2))
