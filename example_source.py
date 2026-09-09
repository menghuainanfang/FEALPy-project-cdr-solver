"""Solve -u''=2 with zero endpoints: f is explicitly provided."""
from cdr_solver import solve_cdr, expression_function as fun
mesh, space, uh, residual, boundary_error = solve_cdr(
    a=fun(1), b=fun(0, vector=True), c=fun(0),
    f=fun(2), gd=fun(0), n=16, p=2)
print('Maximum solution value:', max(uh[:]))
print('Relative residual:', residual)
print('Boundary error:', boundary_error)
