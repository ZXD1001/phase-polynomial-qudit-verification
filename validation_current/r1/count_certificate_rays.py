"""Exact ray inventory for the 23-setting, gap-1/2 full-front-menu certificate.
This is a finite n=4 inventory, not a minimum calibration count or an all-n claim."""
import json
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from lz_full_menu import configs,target_phase
OUT=Path(__file__).resolve().parent
ROOT=((1,0),(0,1),(-1,-1))
def mul(x,y):
 a,b=x;c,d=y;return a*c-b*d,a*d+b*c-b*d
def add(x,y):return x[0]+y[0],x[1]+y[1]
def div(x,y):
 a,b=y;d=a*a-a*b+b*b;z=mul(x,(a-b,-b));return F(z[0],d),F(z[1],d)
def key(v):
 p=next(x for x in v if x!=(0,0));return tuple(div(x,p) for x in v)
paulis=set()
for j in range(3):paulis.add(key([(int(x==j),0) for x in range(3)]))
for t in range(3):
 for j in range(3):paulis.add(key([ROOT[(t*x*x+j*x)%3] for x in range(3)]))
X=list(product(range(3),repeat=4));rows=json.loads((OUT/'lz_fullmenu_exact_half_certificate.json').read_text())['weights']
allrays=set();bysite=[set() for _ in range(4)];example=None
native=set();native_sites=[set() for _ in range(4)]
def qubit_rays(k,s):
 B=[j for j in range(4) if j!=k];out=set()
 for o in product(range(2),repeat=3):
  q=[0,0]
  for x in product(range(2),repeat=4):
   if any(t==0 and x[j]!=oj for j,t,oj in zip(B,s,o)):continue
   ph=sum(x[j]*x[j+1] for j in range(3))-sum(oj*x[j] for j,t,oj in zip(B,s,o) if t)
   q[x[k]]+=(-1)**(ph%2)
  if q==[0,0]:continue
  d=next(x for x in q if x!=0);out.add(tuple(F(x,d) for x in q))
 return out
for row in rows:
 k=row['terminal'];s=row['front_string'];B=[j for j in range(4) if j!=k];qr=qubit_rays(k,s)
 for o in product(range(3),repeat=3):
  q=[(0,0)]*3
  for x in X:
   if any(t==0 and x[j]!=oj for j,t,oj in zip(B,s,o)):continue
   ph=sum(x[j]**2*x[j+1] for j in range(3))-sum((t-1)*x[j]**2+oj*x[j] for j,t,oj in zip(B,s,o) if t)
   q[x[k]]=add(q[x[k]],ROOT[ph%3])
  if all(z==(0,0) for z in q):continue
  ky=key(q);allrays.add(ky);bysite[k].add(ky)
  for kr in qr:native.add((kr,ky));native_sites[k].add((kr,ky))
  if ky not in paulis and example is None:example={'terminal':k,'front_string':s,'qutrit_outcomes':o,'qutrit_conditional_coefficients':q}
res={'n':4,'exact_unique_qutrit_terminal_rays':len(allrays),'non_Pauli_rays':len(allrays-paulis),'per_site_unique_qutrit_rays':[len(s) for s in bysite],'unique_native_six_level_rays':len(native),'native_rays_by_terminal':[len(s) for s in native_sites],'complete_PVM_library_lower_bound_by_terminal':[(len(s)+5)//6 for s in native_sites],'example_non_Pauli_ray':example,'interpretation':'Count of qutrit-factor rays over nonzero branches, Per-terminal PVM lower bounds apply only to exact realization of the specified certificate rank-one rays using complete six-level PVMs and classical coarse graining; not all equivalent protocols.'}
(OUT/'fullmenu_ray_inventory.json').write_text(json.dumps(res,indent=2));print(json.dumps(res,indent=2))
