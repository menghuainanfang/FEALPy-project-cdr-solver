"""1D steady CDR: -(a*u')' + b*u' + c*u = f; source f is required."""
from pathlib import Path
import argparse
import csv
import json
import platform
import warnings
import sys
import importlib.metadata as metadata
import numpy as np
import sympy as sp
from cdr_lfem_solver_1d import CdrLFEMSolver1D
import fealpy
from fealpy.backend import backend_manager as bm
from fealpy.decorator import cartesian
from fealpy.mesh import IntervalMesh
from fealpy.functionspace import LagrangeFESpace
from fealpy.fem import (BilinearForm, LinearForm, ScalarDiffusionIntegrator,
    ScalarConvectionIntegrator, ScalarMassIntegrator, ScalarSourceIntegrator, DirichletBC)

x = sp.Symbol('x', real=True)


def expression_function(expr, vector=False):
    """Convert a SymPy expression to a broadcast-safe FEALPy coefficient."""
    fun = sp.lambdify(x, expr, 'numpy')
    @cartesian
    def evaluate(points):
        xx = points[..., 0]
        values = np.broadcast_to(np.asarray(fun(xx), dtype=float), xx.shape)
        return values[..., None] if vector else values
    return evaluate


def manufactured_case(name):
    if name == 'poisson':
        a, b, c = sp.Integer(1), sp.Integer(0), sp.Integer(0)
    elif name == 'variable':
        a, b, c = 1 + x**2, 1 + x, 2 + x
    elif name == 'reverse':
        a, b, c = 2 + x, 1 - 2*x, 3 + x**2
    elif name == 'cubic':
        a, b, c = 2 + x**3, -1 + x**2, 2 + x
    else:
        raise ValueError(name)
    exact = sp.sin(sp.pi*x) + 1 + x  # Nonzero boundary values: 1 and 2.
    source = -sp.diff(a*sp.diff(exact, x), x) + b*sp.diff(exact, x) + c*exact
    return dict(a=a, b=b, c=c, exact=exact, source=sp.simplify(source))


def _solve(a, b, c, f, gd, n=32, p=1, interval=(0., 1.), q=None, *, formulation):
    """NumPy example adapter around the independently reusable FEALPy core."""
    if bm.get_current_backend().backend_name != 'numpy':
        raise RuntimeError('This SymPy/NumPy example requires numpy; use the core for PyTorch')
    if isinstance(n, bool) or not isinstance(n, int) or n < 1:
        raise ValueError('n must be a positive integer')
    if len(interval) != 2 or not np.all(np.isfinite(interval)) or interval[0] >= interval[1]:
        raise ValueError('Expected finite increasing interval endpoints')
    mesh = IntervalMesh.from_interval_domain(list(interval), nx=n)
    space = LagrangeFESpace(mesh, p=p)
    model = CdrLFEMSolver1D(space, a, b, c, f, gd, q=q, formulation=formulation)
    uh = model.solve()
    return mesh, space, uh, model.relative_residual, model.boundary_error


def solve_transport(a, b, r, f, gd, n=32, p=1, interval=(0., 1.), q=None):
    """Comparison API: (-a*u' + b*u)' + r*u = f; prescribed endpoint values.

    No velocity derivative is required. r denotes the conservative reaction.
    """
    return _solve(a, b, r, f, gd, n, p, interval, q, formulation='conservative')


def solve_cdr(a, b, c, f, gd, n=32, p=1, interval=(0., 1.), q=None):
    """Main API: -(a*u')' + b*u' + c*u = f. f is an independent source input."""
    return _solve(a, b, c, f, gd, n, p, interval, q, formulation='advective')


def equation_balance(mesh, uh, case):
    """Integrated PDE balance using diffusive flux j=-a*u'.

    CG FEM does not guarantee exact cellwise balance. This residual must converge.
    """
    ends = np.array([[0.], [1.]])
    gradients = uh.grad_value(np.array([[1., 0.], [0., 1.]]))
    du = np.array([gradients[0, 0, 0], gradients[-1, 1, 0]])
    values = uh(np.array([[1., 0.], [0., 1.]]))
    uend = np.array([values[0, 0], values[-1, 1]])
    flux = -expression_function(case['a'])(ends)*du
    bcs, weights = mesh.quadrature_formula(8, 'cell').get_quadrature_points_and_weights()
    points = mesh.bc_to_point(bcs)
    reaction = (expression_function(case['b'])(points)*uh.grad_value(bcs)[..., 0]
                + expression_function(case['c'])(points)*uh(bcs))
    source = expression_function(case['source'])(points)
    measure = mesh.entity_measure('cell')
    integral = np.einsum('q,cq,c->', weights, reaction-source, measure)
    abs_integral = np.einsum('q,cq,c->', weights, np.abs(reaction)+np.abs(source), measure)
    imbalance = float(flux[1]-flux[0]+integral)
    scale = max(1., float(np.sum(np.abs(flux))+abs_integral))
    return abs(imbalance), abs(imbalance)/scale


def run(output):
    bm.set_backend('numpy')  # Explicit application-level choice, never on import.
    output.mkdir(parents=True, exist_ok=True)
    rows = []
    last = None
    for name in ('poisson', 'variable', 'reverse', 'cubic'):
        case = manufactured_case(name)
        exact = expression_function(case['exact'])
        gradient = expression_function(sp.diff(case['exact'], x), vector=True)
        for p in (1, 2):
            previous = None
            for n in (8, 16, 32, 64):
                mesh, space, uh, residual, boundary = solve_cdr(
                    expression_function(case['a']), expression_function(case['b'], True),
                    expression_function(case['c']), expression_function(case['source']), exact, n, p)
                l2 = float(mesh.error(exact, uh, q=8))
                h1 = float(mesh.error(gradient, uh.grad_value, q=8))
                balance, relative_balance = equation_balance(mesh, uh, case)
                row = dict(equation_balance=balance, relative_equation_balance=relative_balance, case=name, p=p, n=n, dofs=space.number_of_global_dofs(),
                    L2=l2, H1_seminorm=h1,
                    L2_order=np.log2(previous[0]/l2) if previous else 0.,
                    H1_order=np.log2(previous[1]/h1) if previous else 0.,
                    relative_residual=residual, boundary_error=boundary)
                rows.append(row)
                print(f'{name:8} P{p} n={n:3} L2={l2:.3e} H1semi={h1:.3e} '
                      f'orders={row["L2_order"]:.2f}/{row["H1_order"]:.2f} residual={residual:.2e}')
                assert np.isfinite(l2+h1) and residual < 1e-10 and boundary < 1e-12
                previous = l2, h1
                last = mesh, space, uh, case
            assert row['L2_order'] > p + .8, row
            assert row['H1_order'] > p - .2, row
            assert row['relative_equation_balance'] < .01, row
            data = [r for r in rows if r['case']==name and r['p']==p]
            assert data[-1]['equation_balance'] < data[0]['equation_balance']/4, data
    with (output/'convergence.csv').open('w', newline='', encoding='utf-8') as file:
        writer = csv.DictWriter(file, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    run_metadata = dict(formulation='advective', equation="-(a*u')' + b*u' + c*u = f", reaction_symbol='c', python=platform.python_version(), executable=sys.executable, versions={name: metadata.version(name) for name in ("fealpy", "numpy", "scipy", "sympy", "matplotlib")}, fealpy_path=fealpy.__file__,
                    numpy=np.__version__, cases={k: {s: str(v) for s,v in manufactured_case(k).items()}
                                              for k in ('poisson','variable','reverse','cubic')})
    (output/'run_info.json').write_text(json.dumps(run_metadata, indent=2), encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, grid_axes = plt.subplots(2, 2, figsize=(11, 8))
    axes = grid_axes.ravel()
    mesh, space, uh, case = last
    points = space.interpolation_points()[:, 0]
    idx = np.argsort(points)
    dense = np.linspace(0, 1, 500)
    axes[0].plot(dense, expression_function(case['exact'])(dense[:, None]), label='Exact')
    axes[0].plot(points[idx], uh[:][idx], '.', markersize=3, label='FEM P2, n=64')
    axes[0].set(xlabel='x', ylabel='u', title='Cubic coefficients; u(0)=1, u(1)=2')
    axes[0].legend()
    for p in (1, 2):
        data = [r for r in rows if r['case']=='cubic' and r['p']==p]
        for key in ('L2', 'H1_seminorm'):
            axes[1].loglog([1/r['n'] for r in data], [r[key] for r in data], 'o-', label=f'P{p} {key}')
    axes[1].set(xlabel='h', ylabel='Error', title='Mesh convergence')
    axes[1].legend()
    axes[1].grid(True, which='both', alpha=.3)
    local_t = np.linspace(0, 1, 17)
    bcs = np.stack([1-local_t, local_t], axis=-1)
    samples = mesh.bc_to_point(bcs)[..., 0].ravel()
    numerical = uh(bcs).ravel()
    exact_values = expression_function(case['exact'])(samples[:, None])
    axes[2].plot(samples, numerical-exact_values)
    axes[2].set(xlabel='x', ylabel='u_h - u', title='Pointwise error (Cubic, P2, n=64)')
    for key in ('a', 'b', 'c', 'source'):
        axes[3].plot(dense, expression_function(case[key])(dense[:, None]), label='f (source)' if key=='source' else key)
    axes[3].set(xlabel='x', ylabel='Coefficient', title='Cubic case coefficients')
    axes[3].legend()
    fig.tight_layout()
    fig.savefig(output/'verification.png', dpi=180)
    plt.close(fig)
    print('PASS: convergence, boundary values, residuals and integrated equation balance convergence. Results:', output.resolve())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path(__file__).parent/'results')
    run(parser.parse_args().output)


