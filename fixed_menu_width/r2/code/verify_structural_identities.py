"""Exact finite identities used in the all-size menu-optimality proof.
The all-size conclusion is a mathematical derivation in REPORT.md; these
checks verify its local CRT phases and two-qubit Bell identity, not an
exhaustive search over arbitrary-size protocols.
"""
if not __debug__:
    raise RuntimeError('Run without -O/-OO; verification uses assertions.')
from itertools import product
from pathlib import Path
import json
import sympy as s
BASE=Path(__file__).resolve().parents[1]

def main():
    gate_count=basis_count=0
    for a,b,c,e in product(range(2),range(3),range(2),range(3)):
        x=(3*a+4*b)%6;y=(3*c+4*e)%6
        assert ((2*x*x+3*x)*y-3*a*c-2*b*b*e)%6==0
        gate_count+=1
        for t in range(3):
            assert (2*t*x*x-y*x-3*a*c-2*(t*b*b+b*e))%6==0
            basis_count+=1
    # Unnormalized Bell rays avoid square roots; projectors divide by two.
    rays=[s.Matrix(x) for x in ([1,0,0,1],[1,0,0,-1],[0,1,1,0],[0,1,-1,0])]
    P=[v*v.T/2 for v in rays]
    Y=s.Matrix([[0,-s.I],[s.I,0]]);YY=s.kronecker_product(Y,Y)
    assert P[1]+P[2]-P[0]-P[3]==YY
    def pt0(A):return s.Matrix(4,4,lambda i,j:A[2*(j//2)+i%2,2*(i//2)+j%2])
    assert pt0(YY)==-YY
    for v in ([1,0],[0,1],[1,1],[1,-1]):
        x=s.Matrix(v);q=x*x.T
        assert q.T==q
    out={'status':'EXACT_PASS','CRT_edge_phase_cases':gate_count,'CRT_terminal_basis_phase_cases':basis_count,
      'Bell_identity':'P_Phi- + P_Psi+ - P_Phi+ - P_Psi- = Y tensor Y',
      'partial_transpose_identity':'(Y tensor Y)^T0 = -Y tensor Y',
      'scope':'Exact local identities; all-n theorem additionally uses positivity, perfect completeness, and real-qubit leaf effects as proved in REPORT.md.'}
    (BASE/'results'/'structural_identities_verified.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(out),flush=True)
if __name__=='__main__':main()
