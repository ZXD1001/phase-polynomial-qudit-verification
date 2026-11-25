"""EXACT four-site certificate: canonical Li-Zhu conditional tests using the
same FRONT menu {Z,B0,B1,B2} as the manuscript achieve gap exactly 1/2.
The terminal ray is unrestricted, which is an essential qualification.
Candidate generation may use floating values. --verify-only checks solely
integer/Fraction arithmetic (NumPy object arrays are integer containers).
PSD is checked by an exact rational congruence and strict diagonal dominance.
"""
if not __debug__:
    raise RuntimeError("Exact checks require normal Python; do not use -O or -OO.")

from pathlib import Path
from itertools import product
from fractions import Fraction as F
from math import lcm
import json,sys
import numpy as np
OUT=Path(__file__).resolve().parent;N=81
ROOT=((1,0),(0,1),(-1,-1));X=list(product(range(3),repeat=4))
def add(x,y):return x[0]+y[0],x[1]+y[1]
def mul(x,y):
 a,b=x;c,d=y;return a*c-b*d,a*d+b*c-b*d
def conj(x):return x[0]-x[1],-x[1]
def norm(x):a,b=x;return a*a-a*b+b*b
h=[ROOT[sum(x[i]**2*x[i+1] for i in range(3))%3] for x in X]
def qaccept(k,s,sector):
 B=[j for j in range(4) if j!=k];err=tuple((sector>>(3-j))&1 for j in range(4))
 for a in product(range(2),repeat=4):
  b=tuple(((a[j-1] if j else 0)+(a[j+1] if j<3 else 0))%2 for j in range(4))
  if any((a[j] if t==0 else b[j]) for j,t in zip(B,s)):continue
  if sum(x*y for x,y in zip(a,err))%2:return 0
 return 1

def branch_raw(k,s,o):
 B=[j for j in range(4) if j!=k];q=[(0,0)]*3;phs=[]
 for x,hx in zip(X,h):
  if any(t==0 and x[j]!=oj for j,t,oj in zip(B,s,o)):
   phs.append(None);continue
  ph=sum((t-1)*x[j]**2+oj*x[j] for j,t,oj in zip(B,s,o) if t)%3
  phs.append(ph);q[x[k]]=add(q[x[k]],mul(hx,ROOT[-ph%3]))
 nq=sum(norm(z) for z in q)
 if nq==0:return None
 raw=[(0,0) if ph is None else mul(ROOT[ph],q[x[k]]) for ph,x in zip(phs,X)]
 return np.array([r[0] for r in raw],dtype=object),np.array([r[1] for r in raw],dtype=object),3**sum(t>0 for t in s)*nq

def outer(a,b):
 # (a+b*w)(a+b*w)^dagger
 return np.outer(a,a-b)+np.outer(b,b),-np.outer(a,b)+np.outer(b,a)
def kmat(A,B,C,D):return A@C-B@D,A@D+B@C-B@D

path=OUT/'lz_fullmenu_exact_half_certificate.json'
if '--verify-only' not in sys.argv:
 import sympy as sp
 from lz_full_menu import build
 from scipy.linalg import cholesky,solve_triangular
 d=json.loads((OUT/'lz_fullmenu_half_candidate.json').read_text());ids=[i for i,w in enumerate(d['weights']) if w>1e-10];labs=[d['labels'][i] for i in ids];pp=np.array([d['weights'][i] for i in ids])
 rows=[[1]*len(ids)];rhs=[sp.Rational(1)]
 for sector in range(1,16):
  row=[1-qaccept(k,s,sector) for k,s in labs]
  if abs(float(np.dot(row,pp))-.5)<1e-7:rows.append(row);rhs.append(sp.Rational(1,2))
 A=sp.Matrix(rows);b=sp.Matrix(rhs);aug,piv=A.row_join(b).rref();piv=[j for j in piv if j<len(ids)];free=[j for j in range(len(ids)) if j not in piv]
 p=[sp.Rational(round(float(w)*10**10),10**10) for w in pp]
 for rr,j in enumerate(piv):p[j]=aug[rr,-1]-sum(aug[rr,k]*p[k] for k in free)
 assert A*sp.Matrix(p)==b and all(v>0 for v in p)
 cf={'n':4,'class':'All front strings over {Z,B0,B1,B2}, arbitrary terminal ray, arbitrary terminal/setting weights','proposed_gap':'1/2','weights':[{'terminal':k,'front_string':s,'weight':str(v)} for (k,s),v in zip(labs,p)],'generation_note':'Float weights rationalized under exact qubit active constraints; float inverse Cholesky rounded only to propose congruence.'}
 path.write_text(json.dumps(cf,indent=2))
c=json.loads(path.read_text());settings=[(r['terminal'],r['front_string'],F(r['weight'])) for r in c['weights']]
assert all(p>0 for k,s,p in settings) and sum(p for k,s,p in settings)==1
qrej=[sum(p*(1-qaccept(k,s,j)) for k,s,p in settings) for j in range(1,16)]
assert min(qrej)==F(1,2)
# Build exact M = R3 - I/2 + (3/2)|H><H|.
branches=[];den=54
for k,s,p in settings:
 for o in product(range(3),repeat=3):
  raw=branch_raw(k,s,o)
  if raw is None:continue
  aa,bb,nq=raw;w=p/nq;den=lcm(den,w.denominator);branches.append((aa,bb,w))
AA=np.zeros((N,N),dtype=object);BB=np.zeros((N,N),dtype=object)
for a,b,w in branches:
 oa,ob=outer(a,b);fac=w.numerator*(den//w.denominator);AA-=fac*oa;BB-=fac*ob
for j in range(N):AA[j,j]+=den//2
ha=np.array([z[0] for z in h],dtype=object);hb=np.array([z[1] for z in h],dtype=object);pa,pb=outer(ha,hb);AA+=(den//54)*pa;BB+=(den//54)*pb
assert np.array_equal(AA.T,AA-BB) and np.array_equal(BB.T,-BB)
xa=AA@ha-BB@hb;xb=AA@hb+BB@ha-BB@hb
assert np.array_equal(xa,den*ha) and np.array_equal(xb,den*hb)
if '--verify-only' not in sys.argv:
 from scipy.linalg import cholesky,solve_triangular,eigh
 z=(np.asarray(AA,dtype=float)+np.asarray(BB,dtype=float)*np.exp(2j*np.pi/3))/den
 L=cholesky(z,lower=True);K=solve_triangular(L.conj().T,np.eye(N),lower=False);scale=10**8
 Bk=np.rint(K.imag*2/np.sqrt(3)*scale).astype(np.int64);Ak=np.rint((K.real+K.imag/np.sqrt(3))*scale).astype(np.int64)
 c['congruence']={'scale':scale,'a':Ak.tolist(),'b':Bk.tolist()};c['floating_min_M_eigenvalue_proposal_only']=float(eigh(z,eigvals_only=True)[0]);path.write_text(json.dumps(c,indent=2))
ck=c['congruence'];KA=np.array(ck['a'],dtype=object);KB=np.array(ck['b'],dtype=object)
C,D=kmat(AA,BB,KA,KB);EA,EB=kmat((KA-KB).T,-KB.T,C,D)
assert np.array_equal(EA.T,EA-EB) and np.array_equal(EB.T,-EB)
# |a+b*w| <= |a|+|b|. Strict diagonal dominance implies positive definite.
margins=[]
for j in range(N):
 assert EB[j,j]==0
 margin=int(EA[j,j])-sum(abs(int(EA[j,k]))+abs(int(EB[j,k])) for k in range(N) if k!=j)
 assert margin>0,(j,margin)
 margins.append(margin)
md=F(min(margins),den*int(ck['scale'])**2)
result={'status':'EXACT_PASS','gap':'1/2','number_nonzero_settings':len(settings),'front_menu':'Z,B0,B1,B2','terminal':'canonical conditional ray, NOT restricted to the three terminal bases','qubit_rejection_rates':[str(x) for x in qrej],'M_definition':'R3-I/2+(3/2)P_H4','target_check':'M|H4>=|H4> exactly','PSD_certificate':'K^dagger M K strictly diagonally dominant with positive diagonal, all arithmetic integer; hence M positive definite. Orthogonal qutrit rejection strictly exceeds 1/2.','minimum_scaled_Gershgorin_margin':str(md),'margin_decimal_display':float(md),'exact_gap_upper_reason':'Two endpoint-adjacent qubit error sectors give class-wide 1/2 upper bound, proved separately.'}
(OUT/'lz_fullmenu_exact_half_result.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
