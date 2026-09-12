"""Continuous Lagrange FEM for steady and transient CDR equations."""
from pathlib import Path
from copy import copy
import numpy as np
from fealpy.backend import backend_manager as bm
from fealpy.mesh import IntervalMesh, TriangleMesh, TetrahedronMesh
from fealpy.fem import (BilinearForm, LinearForm, ScalarDiffusionIntegrator,
                       ScalarConvectionIntegrator, ScalarMassIntegrator,
                       ScalarSourceIntegrator, DirichletBC)
from fealpy.solver import spsolve


def box_mesh(pde, n):
    if pde.dim == 1:
        return IntervalMesh.from_interval_domain(pde.box, nx=n)
    if pde.dim == 2:
        return TriangleMesh.from_box(pde.box, nx=n, ny=n)
    if pde.dim == 3:
        return TetrahedronMesh.from_box(pde.box, nx=n, ny=n, nz=n)
    raise ValueError('Mesh construction supports dimensions 1, 2 and 3')


def direct_solve(A, F):
    return spsolve(A, F, solver='scipy')


class CDRLFEMSolver:
    """NumPy/CPU solver with full Dirichlet boundaries on a box."""
    def __init__(self, space, pde, *, q=None, linear_solver=None):
        if bm.get_current_backend().backend_name != 'numpy':
            raise ValueError('This solver uses the NumPy backend')
        self.space, self.mesh, self.pde = space, space.mesh, pde
        self.q = space.p + 4 if q is None else q
        self.linear_solver = direct_solve if linear_solver is None else linear_solver
        points = space.interpolation_points()
        self.boundary = np.zeros(points.shape[0], dtype=bool)
        for i in range(pde.dim):
            self.boundary |= np.isclose(points[:, i], pde.box[2*i], rtol=0, atol=1e-12)
            self.boundary |= np.isclose(points[:, i], pde.box[2*i+1], rtol=0, atol=1e-12)
        self.history = []

    def _check_backend(self):
        if bm.get_current_backend().backend_name != 'numpy':
            raise RuntimeError('Select NumPy before reusing this NumPy solver')

    def check_coefficients(self, t, transient=False):
        self._check_backend()
        bcs, _ = self.mesh.quadrature_formula(self.q, 'cell').get_quadrature_points_and_weights()
        points = self.mesh.bc_to_point(bcs)
        if np.any(self.pde.diffusion(points, t) <= 0):
            raise ValueError('Diffusion must be positive at quadrature points')
        if transient and np.any(self.pde.capacity(points, t) <= 0):
            raise ValueError('Time capacity must be positive at quadrature points')

    def init_integrators(self, t):
        at = lambda f: self.pde.at_time(f, t)
        return (ScalarDiffusionIntegrator(coef=at(self.pde.diffusion), q=self.q),
                ScalarConvectionIntegrator(coef=at(self.pde.convection), q=self.q),
                ScalarMassIntegrator(coef=at(self.pde.reaction), q=self.q),
                ScalarSourceIntegrator(source=at(self.pde.source), q=self.q))

    def linear_system(self, t=0):
        self._check_backend()
        diffusion, convection, reaction, source = self.init_integrators(t)
        form = BilinearForm(self.space)
        form.add_integrator(diffusion, convection, reaction)
        rhs = LinearForm(self.space)
        rhs.add_integrator(source)
        return form.assembly(), rhs.assembly()

    def mass_matrix(self, t):
        self._check_backend()
        form = BilinearForm(self.space)
        form.add_integrator(ScalarMassIntegrator(coef=self.pde.at_time(self.pde.capacity, t), q=self.q))
        return form.assembly()

    def _solve_system(self, A, F, t):
        self.A, self.F = DirichletBC(self.space, gd=self.pde.at_time(self.pde.dirichlet, t),
                                   threshold=self.boundary).apply(A, F)
        uh = self.space.function()
        uh[:] = self.linear_solver(self.A.copy(), self.F.copy())
        self.uh = uh
        residual = np.linalg.norm(self.A @ uh[:] - self.F)/max(1, np.linalg.norm(self.F))
        points = self.space.interpolation_points()[self.boundary]
        boundary = np.max(np.abs(uh[:][self.boundary] - self.pde.dirichlet(points, t)))
        self.history.append(dict(time=float(t), residual=float(residual), boundary=float(boundary)))
        return uh

    def solve(self):
        self.check_coefficients(0)
        return self._solve_system(*self.linear_system(), 0)

    def solve_time(self, end_time, steps, output=None):
        """Backward Euler for d*u_t; initial data are specified at t=0."""
        self._check_backend()
        if end_time <= 0 or steps < 1:
            raise ValueError('Positive final time and step count are required')
        if self.pde.initial is None:
            raise ValueError('Transient solves require initial data')
        uh = self.space.interpolate(self.pde.initial)
        dt = end_time/steps
        entries = []
        if output is not None:
            output = Path(output)
            output.mkdir(parents=True, exist_ok=True)
            self.export_vtu(uh, output/'solution_0000.vtu')
            entries.append((0., 'solution_0000.vtu'))
        for k in range(1, steps+1):
            t = k*dt
            self.check_coefficients(t, transient=True)
            A, F = self.linear_system(t)
            M = self.mass_matrix(t)
            uh = self._solve_system(M + dt*A, M @ uh[:] + dt*F, t)
            if output is not None:
                name = f'solution_{k:04d}.vtu'
                self.export_vtu(uh, output/name)
                entries.append((t, name))
        if output is not None:
            lines = ['<?xml version="1.0"?>', '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian"><Collection>']
            lines += [f'<DataSet timestep="{t:.16g}" group="" part="0" file="{name}"/>' for t, name in entries]
            lines += ['</Collection></VTKFile>']
            (output/'solution.pvd').write_text('\n'.join(lines), encoding='utf-8')
        return uh

    def export_vtu(self, uh, filename):
        """Export vertex samples; high-order interior variation is not exported."""
        self._check_backend()
        bcs = np.eye(self.pde.dim+1)
        samples = np.asarray(uh(bcs))
        nodal = np.empty(self.mesh.number_of_nodes())
        nodal[np.asarray(self.mesh.cell)] = samples
        exported_mesh = copy(self.mesh)
        exported_mesh.nodedata = dict(self.mesh.nodedata, u=nodal)
        exported_mesh.to_vtk(fname=str(filename))
