"""Solve -u''=2 on [0,1] with homogeneous Dirichlet boundaries."""
import numpy as np
from fealpy.backend import backend_manager as bm
from cdr import SymbolicCDR, CDRSolver, box_mesh

bm.set_backend("numpy")
pde = SymbolicCDR(dim=1, exact=["x*(1-x)"])
mesh = box_mesh(pde, n=2)
solver = CDRSolver(pde, mesh, p=1, stabilization="galerkin")
uh = solver.solve_steady()
points = solver.scalar_space.interpolation_points()
values = np.asarray(uh[:]).reshape(1, -1)[0]
order = np.argsort(points[:, 0])
print("x =", points[order, 0])
print("u =", values[order])
np.testing.assert_allclose(values, points[:, 0] * (1-points[:, 0]), atol=1e-12)
print("PASS: three nodal values agree with the hand calculation")
