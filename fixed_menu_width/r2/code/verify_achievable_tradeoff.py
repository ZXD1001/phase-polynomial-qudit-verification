"""Check an analytic achievable width/gap tradeoff.

Exact sector enumeration uses Fraction. Independent n=4 native-six-level
branch construction uses exact cyclotomic support and floating matrices;
its residual and spectrum are cross-checks, not the all-n proof.
"""
if not __debug__:
    raise RuntimeError("Run without -O/-OO.")
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
from itertools import product
from fractions import Fraction as Q
from math import comb
import json,time
import numpy as np
from scipy.linalg import eigvalsh
ROOT6=np.array([[1,0],[0,1],[-1,1],[-1,0],[0,-1],[1,-1]],dtype=np.int64)
BASE=Path(__file__).resolve().parents[1]

def choose(n,k): return comb(n,k) if 0<=k<=n else 0

def exact_sector_check(n,r):
    m=[sum(v%2==c for v in range(n-1)) for c in range(2)]
    d=[min(r,x) for x in m]
    qc=[Q(1) if not x else Q(1,3)+Q(2,3)*y/x for x,y in zip(m,d)]
    weights=[qc[1]/sum(qc),qc[0]/sum(qc)]
    best=Q(1);witness=None
    last=(n-1)%2
    # Binary flags say whether any qubit character is nonzero in each color.
    for s0,s1,z0,z1,zlast in product(range(m[0]+1),range(m[1]+1),range(2),range(2),range(2)):
        if not(s0 or s1 or z0 or z1 or zlast):continue
        ss=[s0,s1];zz=[z0,z1];passes=[]
        for c in range(2):
            passes.append(Q(0) if zz[c] or (c==last and zlast) else Q(choose(m[c]-ss[c],d[c]),choose(m[c],d[c]))*Q(2,3)**ss[c])
        rej=1-sum(p*t for p,t in zip(weights,passes))
        if rej<best:best=rej;witness=[s0,s1,z0,z1,zlast]
    formula=qc[0]*qc[1]/sum(qc)
    assert best==formula,(n,r,best,formula)
    return {'n':n,'width':r,'nonlast_color_counts':m,'adaptive_counts':d,'single_error_color_rejections':list(map(str,qc)),
      'color_weights':list(map(str,weights)),'gap':str(best),'saturating_sector_flags':witness}

def kronall(xs):
    ans=np.ones((1,1),complex)
    for x in xs:ans=np.kron(ans,x)
    return ans

def phase_and_support(n,s):
    X=np.array(list(product(range(6),repeat=n)),int)
    ph=(((2*X[:,:-1]**2+3*X[:,:-1])*X[:,1:]).sum(axis=1))%6
    aa=ROOT6[ph,0].reshape((6,)*n);bb=ROOT6[ph,1].reshape((6,)*n)
    for k,code in enumerate(s):
        if code==0:continue
        a=np.moveaxis(aa,k,0);b=np.moveaxis(bb,k,0);oa=np.zeros_like(a);ob=np.zeros_like(b)
        for j in range(6):
            for x in range(6):
                c,d=ROOT6[(-2*(code-1)*x*x+j*x)%6]
                oa[j]+=c*a[x]-d*b[x]
                ob[j]+=d*a[x]+(c+d)*b[x]
        aa=np.moveaxis(oa,0,k);bb=np.moveaxis(ob,0,k)
    nums=(aa*aa+aa*bb+bb*bb).ravel()
    assert int(nums.sum())==6**(n+sum(t!=0 for t in s))
    return nums!=0

def native_n4():
    tic=time.time();n=4;D=6**n;X=np.array(list(product(range(6),repeat=n)),int)
    w=np.exp(2j*np.pi/6);x=np.arange(6);F=w**(-np.outer(x,x))/np.sqrt(6)
    bases=[np.eye(6,dtype=complex)]+[w**(2*t*x*x)[:,None]*F for t in range(3)]
    # seven top-level templates: odd block 2/5, six partially adaptive even blocks 1/10.
    templates=[(1,{0:0,2:0,3:1},Q(2,5))]
    for t in range(3):
        templates.extend([(0,{1:0,2:1+t,3:0},Q(1,10)),(2,{0:1+t,1:0,3:0},Q(1,10))])
    Omega=np.zeros((D,D),complex);max_orth_res=0.;ranks=[]
    for terminal,front,weight in templates:
        cols=[]
        for t in range(3):
            s=[front.get(v,1+t) for v in range(n)]
            U=kronall([bases[b] for b in s])
            selection=(X[:,terminal+1]%3==t)&phase_and_support(n,s)
            cols.append(U[:,selection])
        V=np.concatenate(cols,axis=1);T=V@V.conj().T
        max_orth_res=max(max_orth_res,float(np.linalg.norm(V.conj().T@V-np.eye(V.shape[1]),ord='fro')))
        ranks.append(V.shape[1]);Omega+=float(weight)*T
    ph=(((2*X[:,:-1]**2+3*X[:,:-1])*X[:,1:]).sum(axis=1))%6;diag=w**ph
    S=np.roll(np.eye(6),1,axis=0);A=(np.eye(6)+np.linalg.matrix_power(S,3))/2
    C=(np.eye(6)+np.linalg.matrix_power(S,2)+np.linalg.matrix_power(S,4))/3;P=np.ones((6,6))/6;I=np.eye(6)
    expected=.4*kronall([I,P,I,P])+.3*kronall([P,I,A@(2*I+C)/3,I])+.3*kronall([A@(2*I+C)/3,I,P,I])
    expected=diag[:,None]*expected*diag.conj()[None,:]
    residual=float(np.linalg.norm(Omega-expected,ord='fro'))
    psi=diag/np.sqrt(D);accept=float(np.linalg.norm(Omega@psi-psi))
    vals=eigvalsh(Omega,check_finite=False);gap=1-float(vals[-2])
    assert residual<1e-10 and abs(gap-.4)<1e-10 and accept<1e-10
    return {'status':'NATIVE_FLOATING_CROSSCHECK_PASS_WITH_EXACT_SUPPORT','n':4,'templates':7,'template_ranks':ranks,
      'operator_formula_Frobenius_residual':residual,'target_residual':accept,'accepted_columns_orthogonality_residual':max_orth_res,
      'gap_float':gap,'exact_gap_from_analytic_sector_proof':'2/5','seconds':time.time()-tic}

if __name__=='__main__':
    checks=[exact_sector_check(n,r) for n in range(2,21) for r in range((n+1)//2+1)]
    native=native_n4()
    result={'status':'PASS','analytic_formula':'q_c = 1/3+(2/3)*min(r,m_c)/m_c (or 1 for m_c=0); nu=q0*q1/(q0+q1)',
      'exact_sector_checks':checks,'native_n4':native,'evidence_scope':'Sector calculation verifies the algebraic construction; proof for all n is in REPORT.md. Floating native n4 matrix residual is a cross-check only.'}
    (BASE/'results'/'achievable_tradeoff_verified.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(native),flush=True)
    print('EXACT_SECTOR_CHECKS',len(checks),flush=True)
