"""Independent reconstruction of SI (S97)--(S104), not the authors' software.
Exact arithmetic: matrices A+B*w over Q, w^2+w+1=0. No floating support tests.
Run: python exact_certificate.py
"""
if not __debug__:
    raise RuntimeError("Exact checks require normal Python; do not use -O or -OO.")

from itertools import product
from pathlib import Path
import json
import sympy as s

class K:
    def __init__(self,A,B=None):
        self.A=s.Matrix(A); self.B=s.zeros(*self.A.shape) if B is None else s.Matrix(B)
    def __add__(self,o): return K(self.A+o.A,self.B+o.B)
    def __sub__(self,o): return K(self.A-o.A,self.B-o.B)
    def __matmul__(self,o):
        return K(self.A*o.A-self.B*o.B,self.A*o.B+self.B*o.A-self.B*o.B)
    def scale(self,c): return K(c*self.A,c*self.B)
    def zero(self): return self.A==s.zeros(*self.A.shape) and self.B==s.zeros(*self.B.shape)
    def tr(self): return (s.trace(self.A),s.trace(self.B))
    def dump(self): return {'a':[[str(v) for v in row] for row in self.A.tolist()], 'b':[[str(v) for v in row] for row in self.B.tolist()]}
ROOT=[(1,0),(0,1),(-1,-1)]
I=K(s.eye(9)); ZERO=K(s.zeros(9))
configs=list(product(range(3),repeat=2))

def phase_col(setting,outcome,x):
    z=0
    for a,j,y in zip(setting,outcome,x):
        if a==0:
            if j!=y: return None
        else: z+=(a-1)*y*y+j*y
    return z%3

def projector(rawphases,denom):
    A=s.zeros(9);B=s.zeros(9)
    for i,ph in enumerate(rawphases):
        if ph is None: continue
        for j,qh in enumerate(rawphases):
            if qh is None: continue
            a,b=ROOT[(ph-qh)%3]
            A[i,j]=s.Rational(a,denom);B[i,j]=s.Rational(b,denom)
    return K(A,B)

reject={}; supports={}
for setting in product(range(4),repeat=2):
    T=ZERO; good=[]
    for out in product(range(3),repeat=2):
        col=[phase_col(setting,out,x) for x in configs]
        coeff=[0,0]
        for (x,y),ph in zip(configs,col):
            if ph is None:continue
            a,b=ROOT[(x*x*y-ph)%3];coeff[0]+=a;coeff[1]+=b
        if coeff==[0,0]:continue
        good.append(out)
        T=T+projector(col,3**sum(a>0 for a in setting))
    assert ((T@T)-T).zero()
    reject[setting]=I-T;supports[str(setting)]=good
P=projector([(x*x*y)%3 for x,y in configs],9)
p={(0,1):s.Rational(2,11)}
p.update({(a,b):s.Rational(1,11) for a in (1,2,3) for b in (0,2,3)})
R=ZERO
for k,v in p.items():R=R+reject[k].scale(v)
poly=R
for k in (3,4,5):poly=poly@(R-I.scale(s.Rational(k,11)))
assert poly.zero()
ker=I
for k in (3,4,5):ker=ker@(I-R.scale(s.Rational(11,k)))
assert (ker-P).zero()
E=((R-I.scale(s.Rational(4,11)))@(R-I.scale(s.Rational(5,11)))).scale(s.Rational(121,2))-P.scale(10)
L=reject[(0,1)];EL=E@L
assert ((E@E)-E).zero() and E.tr()==(4,0)
assert ((E@L)-(L@E)).zero()
assert ((EL@EL)-EL).zero() and EL.tr()==(2,0)
sigma=(E-EL).scale(s.Rational(4,11))+EL.scale(s.Rational(3,22))
assert sigma.tr()==(1,0) and (P@sigma).zero()
expect={}
for k,Q in reject.items():
    v=(Q@sigma).tr()
    target=s.Rational(3,11) if k in p else s.Rational(3,22) if k[0] in (1,2,3) and k[1]==1 else s.Rational(0)
    assert v==(target,0),(k,v,target)
    expect[str(k)]=str(v[0])
# Spectral projectors independently give multiplicities of all eigenvalues.
mult={}
for a in (0,3,4,5):
    Q=I
    for b in (0,3,4,5):
        if a!=b: Q=Q@((R-I.scale(s.Rational(b,11))).scale(s.Rational(11,a-b)))
    assert ((Q@Q)-Q).zero()
    mult[str(s.Rational(a,11))]=str(Q.tr()[0])
res={'status':'EXACT_PASS','arithmetic':'Q(w3), w3^2+w3+1=0; SymPy rational matrices',
     'original_software_used':False,'input':'H2 amplitude w3^(x^2*y)/3; local bases SI S97',
     'checks':['all 16 support PVM projectors idempotent','S100','S101','E projector trace4','EL projector trace2 and [E,L]=0','sigma positive by orthogonal projector decomposition, trace1, target orthogonal','all16 S103 expectations'],
     'weights':{str(k):str(v) for k,v in p.items()},'supports':supports,
     'R_eigenvalue_multiplicities':mult,'dual_expectations':expect,'R':R.dump(),'sigma':sigma.dump()}
out=Path(__file__).with_name('exact_certificate_result.json');out.write_text(json.dumps(res,indent=2),encoding='utf-8')
print(json.dumps({k:res[k] for k in ('status','R_eigenvalue_multiplicities','dual_expectations')},indent=2))
