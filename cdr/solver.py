"""Block Lagrange FEM, SUPG and BDF2 for vector CDR systems."""

from pathlib import Path
from copy import copy
from fealpy.backend import backend_manager as bm
from fealpy.mesh import IntervalMesh, TriangleMesh, TetrahedronMesh
from fealpy.functionspace import LagrangeFESpace, TensorFunctionSpace
from fealpy.fem import (
    BilinearForm,
    LinearForm,
    BlockForm,
    DirichletBC,
    ScalarDiffusionIntegrator,
    ScalarConvectionIntegrator,
    ScalarMassIntegrator,
    ScalarSourceIntegrator,
)
from fealpy.fem.integrator import ConstIntegrator
from fealpy.solver import spsolve
from .pde import at_time
from .stabilization import supg_terms


def box_mesh(pde, n):
    if pde.dim == 1:
        return IntervalMesh.from_interval_domain(pde.domain(), nx=n)
    if pde.dim == 2:
        return TriangleMesh.from_box(pde.domain(), nx=n, ny=n)
    return TetrahedronMesh.from_box(pde.domain(), nx=n, ny=n, nz=n)


class CDRSolver:
    """Scalar (m=1) or coupled vector solver; NumPy/CPU, affine simplices P1/P2."""

    def __init__(self, pde, mesh, p=1, *, stabilization="supg", q=None, linear_solver=None):
        self._check_backend()
        if p not in (1, 2) or stabilization not in ("supg", "galerkin"):
            raise ValueError("Use P1/P2 and supg/galerkin")
        if mesh.geo_dimension() != pde.dim or mesh.top_dimension() != pde.dim:
            raise ValueError("PDE and full-dimensional mesh must agree")
        self.pde, self.mesh = pde, mesh
        self.scalar_space = LagrangeFESpace(mesh, p=p)
        self.space = TensorFunctionSpace(self.scalar_space, shape=(pde.components, -1))
        self.stabilization, self.q = stabilization, p + 4 if q is None else q
        self.linear_solver = linear_solver
        scalar_boundary = pde.is_dirichlet_boundary(self.scalar_space.interpolation_points())
        self.boundary = bm.tile(scalar_boundary, pde.components)
        self.history = []

    @staticmethod
    def _check_backend():
        if bm.get_current_backend().backend_name != "numpy":
            raise RuntimeError("The validated solver and sparse direct engine require NumPy/CPU")

    def check_coefficients(self, t):
        self._check_backend()
        bcs, _ = self.mesh.quadrature_formula(self.q, "cell").get_quadrature_points_and_weights()
        p = self.mesh.bc_to_point(bcs)
        for field in (self.pde.diffusion_coef, self.pde.capacity_coef):
            value = field(p, t)
            if not bool(bm.all(bm.isfinite(value))) or bool(bm.any(value <= 0)):
                raise ValueError(
                    "Diffusion and capacity must be finite and positive at quadrature points"
                )

    def assemble(self, t):
        """Assemble all coefficients at the requested time; no stale time cache."""
        self.check_coefficients(t)
        space, m, q = self.scalar_space, self.pde.components, self.q
        dofs = space.cell_to_dof()
        if self.stabilization == "supg":
            blocks, mass, loads = supg_terms(space, self.pde, t, q)
        forms = []
        rhs = []
        massforms = [[None for _ in range(m)] for _ in range(m)]
        for i in range(m):
            row = []
            for j in range(m):
                form = BilinearForm(space)
                form.add_integrator(
                    ScalarMassIntegrator(coef=at_time(self.pde.reaction_coef, t, (i, j)), q=q)
                )
                if i == j:
                    form.add_integrator(
                        ScalarDiffusionIntegrator(coef=at_time(self.pde.diffusion_coef, t), q=q)
                    )
                    form.add_integrator(
                        ScalarConvectionIntegrator(coef=at_time(self.pde.convection_coef, t), q=q)
                    )
                if self.stabilization == "supg":
                    form.add_integrator(ConstIntegrator(blocks[i][j], dofs))
                row.append(form)
            forms.append(row)
            load = LinearForm(space)
            load.add_integrator(
                ScalarSourceIntegrator(source=at_time(self.pde.source, t, (i,)), q=q)
            )
            mf = BilinearForm(space)
            mf.add_integrator(ScalarMassIntegrator(coef=at_time(self.pde.capacity_coef, t), q=q))
            if self.stabilization == "supg":
                load.add_integrator(ConstIntegrator(loads[i], dofs))
                mf.add_integrator(ConstIntegrator(mass, dofs))
            rhs.append(load.assembly())
            massforms[i][i] = mf
        return BlockForm(forms).assembly(), BlockForm(massforms).assembly(), bm.concatenate(rhs)

    def _solve(self, A, F, t):
        self.A, self.F = DirichletBC(
            self.space, gd=at_time(self.pde.dirichlet, t), threshold=self.boundary
        ).apply(A, F)
        if self.linear_solver is None:
            values = spsolve(self.A.copy(), self.F.copy(), solver="scipy")
        else:
            values = self.linear_solver(self.A.copy(), self.F.copy())
        if values.shape != self.F.shape or not bool(bm.all(bm.isfinite(values))):
            raise ValueError("Linear solver must return a finite vector with the system shape")
        self.uh = self.space.function()
        self.uh[:] = values
        residual = float(bm.linalg.norm(self.A @ values - self.F)) / max(
            1.0, float(bm.linalg.norm(self.F))
        )
        exact_boundary = self.pde.dirichlet(self.scalar_space.interpolation_points(), t).T.reshape(
            -1
        )
        boundary_error = float(
            bm.max(bm.abs(values[self.boundary] - exact_boundary[self.boundary]))
        )
        self.history.append(dict(time=float(t), residual=residual, boundary_error=boundary_error))
        return self.uh

    def solve_steady(self, t=0.0):
        self.history = []
        A, _, F = self.assemble(t)
        return self._solve(A, F, t)

    def solve_time(self, end_time, steps, *, scheme="bdf2", output=None):
        """Uniform-step BDF2 with one backward-Euler start; optional BE comparison."""
        self._check_backend()
        if (
            end_time <= 0
            or not isinstance(steps, int)
            or steps < 1
            or scheme not in ("bdf2", "backward_euler")
        ):
            raise ValueError("Use positive end_time, integer steps and bdf2/backward_euler")
        self.history = []
        points = self.scalar_space.interpolation_points()
        old = self.pde.initial(points).T.reshape(-1).copy()
        initial_boundary = self.pde.dirichlet(points, 0.0).T.reshape(-1)
        old[self.boundary] = initial_boundary[self.boundary]
        previous = None
        dt = end_time / steps
        frames = []
        if output is not None:
            output = Path(output)
            output.mkdir(parents=True, exist_ok=True)
            first = self.space.function()
            first[:] = old
            self.export_vtu(first, output / "solution_0000.vtu")
            frames.append((0.0, "solution_0000.vtu"))
        for k in range(1, steps + 1):
            t = k * dt
            A, M, F = self.assemble(t)
            if scheme == "bdf2" and k > 1:
                left, right = 1.5 * M + dt * A, M @ (2 * old - 0.5 * previous) + dt * F
            else:
                left, right = M + dt * A, M @ old + dt * F
            uh = self._solve(left, right, t)
            previous, old = old, uh[:].copy()
            if output is not None:
                name = f"solution_{k:04d}.vtu"
                self.export_vtu(uh, output / name)
                frames.append((t, name))
        if output is not None:
            lines = [
                '<?xml version="1.0"?>',
                '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian"><Collection>',
            ]
            lines += [
                f'<DataSet timestep="{t:.16g}" group="" part="0" file="{name}"/>'
                for t, name in frames
            ]
            lines += ["</Collection></VTKFile>"]
            (output / "solution.pvd").write_text("\n".join(lines), encoding="utf-8")
        return uh

    def errors(self, t=0.0):
        l2 = self.mesh.error(at_time(self.pde.solution, t), self.uh, q=self.q + 2)
        h1 = self.mesh.error(at_time(self.pde.gradient, t), self.uh.grad_value, q=self.q + 2)
        return float(l2), float(h1)

    def export_vtu(self, uh, filename):
        """Vertex samples; individual components and a VTK vector for m<=3."""
        self._check_backend()
        values = bm.zeros((self.mesh.number_of_nodes(), self.pde.components), dtype=self.mesh.ftype)
        values[self.mesh.cell] = uh(bm.eye(self.pde.dim + 1))
        mesh = copy(self.mesh)
        mesh.nodedata = dict(self.mesh.nodedata)
        for i in range(self.pde.components):
            mesh.nodedata[f"u_{i}"] = values[:, i]
        mesh.nodedata["u_magnitude"] = bm.linalg.norm(values, axis=-1)
        if self.pde.components <= 3:
            vector = bm.zeros((len(values), 3), dtype=values.dtype)
            vector[:, : self.pde.components] = values
            mesh.nodedata["u"] = vector
        mesh.to_vtk(fname=str(filename))
