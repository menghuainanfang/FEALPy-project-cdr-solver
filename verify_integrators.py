"""Small independent checks of FEALPy claims used by the research audit."""
import numpy as np
from fealpy.mesh import IntervalMesh
from fealpy.functionspace import LagrangeFESpace
from fealpy.fem import ScalarDiffusionIntegrator, ScalarConvectionIntegrator, ScalarMassIntegrator
from cdr_solver import expression_function

mesh = IntervalMesh.from_interval_domain([0, 1], nx=1)
space = LagrangeFESpace(mesh, p=1)
checks = [
    ('diffusion', ScalarDiffusionIntegrator(coef=1, q=4), np.array([[1., -1.], [-1., 1.]])),
    ('convection', ScalarConvectionIntegrator(coef=expression_function(1, True), q=4), np.array([[-.5, .5], [-.5, .5]])),
    ('reaction', ScalarMassIntegrator(coef=1, q=4), np.array([[2., 1.], [1., 2.]])/6),
]
for name, integrator, expected in checks:
    actual = np.asarray(integrator.assembly(space))[0]
    np.testing.assert_allclose(actual, expected, atol=1e-13)
    print(name, 'PASS', actual.tolist())
mesh = IntervalMesh.from_interval_domain([0, 1], nx=8)
for p in (1, 2):
    space = LagrangeFESpace(mesh, p=p)
    points = space.interpolation_points()[:, 0]
    endpoints = np.isclose(points, 0) | np.isclose(points, 1)
    print(f'P{p}: default boundary count={np.count_nonzero(space.is_boundary_dof())}; explicit endpoints={np.count_nonzero(endpoints)}')
