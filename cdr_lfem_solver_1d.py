"""Reusable 1D CDR finite-element solver; no import-time backend changes.

Target integration location: fealpy/fem/cdr_lfem_solver_1d.py.
No model-manager registration or package export is performed by this module.
"""
import warnings
from fealpy.backend import backend_manager as bm
from fealpy.decorator import cartesian
from fealpy.mesh import IntervalMesh
from fealpy.functionspace import LagrangeFESpace
from fealpy.fem import (BilinearForm, LinearForm, ScalarDiffusionIntegrator,
    ScalarConvectionIntegrator, ScalarMassIntegrator, ScalarSourceIntegrator,
    DirichletBC)
from fealpy.solver import spsolve

__all__ = ['CdrLFEMSolver1D']


class CdrLFEMSolver1D:
    """Solve -(a*u')' + b*u' + c*u = f on a connected interval.

    Parameters a/b/c/f/gd are scalars or Cartesian callables. Callables must
    use the active backend. Scalars return (...); b may return (..., 1).
    space is a caller-owned P1/P2 LagrangeFESpace on an IntervalMesh.
    NumPy/CPU supports P1/P2; PyTorch/CPU supports P1 only. The caller selects backend.
    linear_solver, if supplied, has signature (A, F) -> coefficient vector.
    A/F are FEALPy sparse/backend objects. The default uses FEALPy spsolve
    with solver='scipy'. No files, plots, logger configuration or global
    backend changes are produced. The mesh/space and input arrays are not
    modified. Reassemble on each solve; reuse across backend switches is
    explicitly rejected. c is r for the optional conservative formulation.
    """
    def __init__(self, space, a, b, c, f, gd, *, q=None,
                 formulation='advective', linear_solver=None):
        if not isinstance(space, LagrangeFESpace):
            raise TypeError('space must be a LagrangeFESpace')
        if not isinstance(space.mesh, IntervalMesh) or space.mesh.geo_dimension() != 1:
            raise ValueError('Only a one-dimensional IntervalMesh is supported')
        if space.p not in (1, 2):
            raise ValueError('Only P1 and P2 are validated')
        if formulation not in ('advective', 'conservative'):
            raise ValueError('Unknown formulation')
        self.backend = bm.get_current_backend().backend_name
        if self.backend not in ('numpy', 'pytorch'):
            raise ValueError('Validated backends: numpy and pytorch on CPU')
        if self.backend == 'pytorch' and space.p == 2:
            raise NotImplementedError('PyTorch P2 is blocked by IntervalMesh interpolation_points')
        self.space = space
        self.mesh = space.mesh
        points = space.interpolation_points()
        if str(bm.get_device(points)) != 'cpu':
            raise ValueError('The default sparse direct solve is CPU-only')
        if q is not None and (isinstance(q, bool) or not isinstance(q, int) or q < 1):
            raise ValueError('q must be a positive integer')
        self.q = space.p + 4 if q is None else q
        self.formulation = formulation
        self.linear_solver = linear_solver
        self.a, self.b = self._coefficient(a), self._coefficient(b, vector=True)
        self.c, self.f, self.gd = map(self._coefficient, (c, f, gd))
        # Explicit endpoint mask avoids the current IntervalMesh BC defect.
        xx = points[..., 0]
        left, right = bm.min(xx), bm.max(xx)
        if not bool(right > left):
            raise ValueError('The interval must have positive length')
        self.boundary = (xx == left) | (xx == right)
        if int(bm.sum(self.boundary)) != 2:
            raise ValueError('Expected exactly two endpoint degrees of freedom')
        self.A = self.F = self.uh = None

    @staticmethod
    def _coefficient(value, vector=False):
        if callable(value) and getattr(value, 'coordtype', 'cartesian') != 'cartesian':
            raise ValueError('Coefficients must use Cartesian coordinates')
        @cartesian
        def evaluate(points):
            result = value(points) if callable(value) else value
            result = bm.asarray(result, **bm.context(points))
            shape = points.shape[:-1]
            if vector and result.shape == shape + (1,):
                pass
            else:
                result = bm.broadcast_to(result, shape)
                if vector:
                    result = result[..., None]
            if not bool(bm.all(bm.isfinite(result))):
                raise ValueError('Coefficient/source/boundary values must be finite')
            return result
        return evaluate

    def _check_backend(self):
        if bm.get_current_backend().backend_name != self.backend:
            raise RuntimeError('Backend changed after solver construction')

    def linear_system(self):
        """Return a fresh unconstrained FEALPy sparse matrix and load vector."""
        self._check_backend()
        bcs, _ = self.mesh.quadrature_formula(self.q, 'cell').get_quadrature_points_and_weights()
        points = self.mesh.bc_to_point(bcs)
        diffusion, velocity = self.a(points), self.b(points)[..., 0]
        if bool(bm.any(diffusion <= 0)):
            raise ValueError('Diffusion must be positive at quadrature points')
        h = self.mesh.entity_measure('cell')[:, None]
        self.max_peclet = float(bm.max(bm.abs(velocity) * h / (2 * diffusion)))
        if self.max_peclet > 1:
            warnings.warn(f'Mesh Peclet number {self.max_peclet:.3g} > 1: '
                          'unstabilized Galerkin may oscillate', RuntimeWarning)
        form = BilinearForm(self.space)
        form.add_integrator(ScalarDiffusionIntegrator(coef=self.a, q=self.q))
        form.add_integrator(ScalarMassIntegrator(coef=self.c, q=self.q))
        convection = BilinearForm(self.space)
        convection.add_integrator(ScalarConvectionIntegrator(coef=self.b, q=self.q))
        B = convection.assembly()
        A = form.assembly() + (-B.T if self.formulation == 'conservative' else B)
        rhs = LinearForm(self.space)
        rhs.add_integrator(ScalarSourceIntegrator(source=self.f, q=self.q))
        return A, rhs.assembly()

    def apply_bc(self, A, F):
        """Apply prescribed endpoint values using FEALPy DirichletBC."""
        self._check_backend()
        return DirichletBC(self.space, gd=self.gd, threshold=self.boundary).apply(A, F)

    def solve(self):
        """Return a fresh FEALPy Function; retain constrained A/F for inspection."""
        self._check_backend()
        self.A, self.F = self.apply_bc(*self.linear_system())
        # SciPy may sort shared CSR storage in-place; isolate solver inputs.
        solve_A, solve_F = self.A.copy(), bm.copy(self.F)
        values = (spsolve(solve_A, solve_F, solver='scipy') if self.linear_solver is None
                  else self.linear_solver(solve_A, solve_F))
        values = bm.asarray(values, **bm.context(self.F))
        if values.shape != self.F.shape or not bool(bm.all(bm.isfinite(values))):
            raise ValueError('Linear solver returned an invalid coefficient vector')
        uh = self.space.function()
        uh[:] = values
        self.uh = uh
        self.relative_residual = float(bm.linalg.norm(self.A @ values - self.F)) / max(
            1., float(bm.linalg.norm(self.F)))
        points = self.space.interpolation_points()[self.boundary]
        self.boundary_error = float(bm.max(bm.abs(values[self.boundary] - self.gd(points))))
        return uh
