"""Independent checks for -(a*u')' + b*u' + c*u = f."""
import json
import argparse
from pathlib import Path
import warnings
import numpy as np
import sympy as sp
from cdr_solver import x, expression_function as fun, manufactured_case, solve_transport, solve_cdr

output={}
case=manufactured_case('variable')
args=[fun(case['a']),fun(case['b'],True)]
source=fun(case['source']); exact=fun(case['exact'])
comparisons=[]
for p in (1,2):
    main=solve_cdr(*args,fun(case['c']),source,exact,n=32,p=p)
    converted=solve_transport(*args,fun(case['c']-sp.diff(case['b'],x)),source,exact,n=32,p=p)
    error=float(np.max(np.abs(main[2][:]-converted[2][:])))
    assert error<1e-11
    comparisons.append(dict(p=p,max_dof_difference=error))
output['equivalent_forms']=comparisons

# All coefficients active: u=1+x^2, a=1+x^2, b=1+x, c=2+x.
# Hand differentiation gives f=x^3-2*x^2+3*x, independently of manufactured_case.
result=solve_cdr(fun(1+x*x),fun(1+x,True),fun(2+x),fun(x**3-2*x*x+3*x),fun(1+x*x),n=16,p=2)
error=float(np.max(np.abs(result[2][:]-fun(1+x*x)(result[1].interpolation_points()))))
assert error<1e-12
output['independent_polynomial_solution']=dict(p=2,n=16,max_error=error,source='x**3-2*x**2+3*x')

# f must affect the solution: fixed coefficients/boundary, change f=0 to f=2.
zero=solve_cdr(fun(1),fun(0,True),fun(0),fun(0),fun(0),n=16,p=2)
nonzero=solve_cdr(fun(1),fun(0,True),fun(0),fun(2),fun(0),n=16,p=2)
points=nonzero[1].interpolation_points()
error=float(np.max(np.abs(nonzero[2][:]-fun(x*(1-x))(points))))
assert np.max(np.abs(zero[2][:]))<1e-13 and error<1e-12
output['source_f_check']=dict(f0_max=float(np.max(np.abs(zero[2][:]))),f2_max=float(np.max(nonzero[2][:])),f2_exact_error=error)

# Same source but wrongly identifying r with c in the other formulation.
wrong=[]
for n in (16,32,64):
    mesh,space,uh,_,_=solve_transport(*args,fun(case['c']),source,exact,n=n,p=2)
    wrong.append(dict(n=n,L2_error=float(mesh.error(exact,uh,q=8))))
assert wrong[-1]['L2_error']>.01
assert .9<wrong[-1]['L2_error']/wrong[-2]['L2_error']<1.1
output['conservative_without_conversion']=wrong

with warnings.catch_warnings(record=True) as captured:
    warnings.simplefilter('always')
    mesh,space,uh,residual,boundary=solve_cdr(fun(.001),fun(1,True),fun(0),fun(0),fun(x),n=32,p=1)
assert captured and np.min(uh[:])<-.01
output['convection_dominated_limitation']=dict(epsilon=.001,n=32,p=1,Pe=15.625,
    min_nodal_value=float(np.min(uh[:])),max_nodal_value=float(np.max(uh[:])),relative_residual=residual,
    boundary_error=boundary,warning=str(captured[0].message),accepted_as_supported=False)
parser=argparse.ArgumentParser()
parser.add_argument('--output',type=Path,default=Path(__file__).with_name('results'))
output_dir=parser.parse_args().output
output_dir.mkdir(parents=True,exist_ok=True)
output_dir.joinpath('form_verification.json').write_text(json.dumps(output,indent=2),encoding='utf-8')
print(json.dumps(output,indent=2))
print('PASS: independent source, polynomial solution, form conversion and known limitation checks')
