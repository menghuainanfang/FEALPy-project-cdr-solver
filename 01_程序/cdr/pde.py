"""PDE interfaces: d*u_t - div(a*grad(u)) + b.grad(u) + C*u = f."""

import sympy as sp
from fealpy.backend import backend_manager as bm
from fealpy.decorator import cartesian


class CDRModel:
    """Shared positive scalar a,d; spatial velocity b; component reaction C.

    solution/source/dirichlet return (..., components), even for one component.
    gradient returns (..., components, dim). Coefficients depend on (p,t).
    """

    def __init__(self, dim, components=1, box=None):
        if dim not in (1, 2, 3) or components < 1:
            raise ValueError("Use dimensions 1/2/3 and a positive component count")
        self.dim, self.components = dim, components
        self.box = list(box) if box is not None else [0.0, 1.0] * dim

    def domain(self):
        return list(self.box)

    @cartesian
    def is_dirichlet_boundary(self, p):
        mask = bm.zeros(p.shape[:-1], dtype=bm.bool)
        for axis in range(self.dim):
            mask |= bm.abs(p[..., axis] - self.box[2 * axis]) < 1e-12
            mask |= bm.abs(p[..., axis] - self.box[2 * axis + 1]) < 1e-12
        return mask

    @cartesian
    def dirichlet(self, p, t=0.0):
        return self.solution(p, t)

    @cartesian
    def initial(self, p):
        return self.solution(p, 0.0)


def at_time(function, t, component=None):
    """Bind time while preserving the Cartesian coefficient convention."""

    @cartesian
    def evaluate(p):
        value = function(p, t)
        return value if component is None else value[(...,) + component]

    return evaluate


class SymbolicCDR(CDRModel):
    """Optional NumPy symbolic data; public PDE functions remain class methods.

    Use x,y,z,t symbols. reaction is an m-by-m matrix; scalar reaction means c*I.
    Supply source/boundary/initial for physical data without an exact solution.
    """

    def __init__(
        self,
        dim,
        *,
        exact=None,
        source=None,
        boundary=None,
        initial=None,
        a=1,
        b=None,
        reaction=0,
        d=1,
        components=None,
        box=None,
    ):
        count = len(exact) if exact is not None else components
        if count is None:
            raise ValueError("Specify components when no exact solution is supplied")
        super().__init__(dim, count, box)
        x, y, z, t = sp.symbols("x y z t")
        coordinates = (x, y, z)[:dim]
        symbols = {str(v): v for v in (x, y, z, t)}

        def expr(value):
            value = sp.sympify(value)
            return value.xreplace(
                {v: symbols[str(v)] for v in value.free_symbols if str(v) in symbols}
            )

        a, d = expr(a), expr(d)
        b = [expr(v) for v in (b if b is not None else [0] * dim)]
        C = (
            sp.Matrix(reaction)
            if isinstance(reaction, (list, tuple, sp.MatrixBase))
            else sp.eye(count) * expr(reaction)
        )
        C = C.applyfunc(expr)
        if C.shape != (count, count) or len(b) != dim:
            raise ValueError("Coefficient dimensions do not match the PDE")
        u = [expr(v) for v in exact] if exact is not None else None
        if source is None:
            if u is None:
                raise ValueError("Provide source or exact solution")
            source = [
                d * sp.diff(u[i], t)
                + sum(
                    -sp.diff(a * sp.diff(u[i], v), v) + b[j] * sp.diff(u[i], v)
                    for j, v in enumerate(coordinates)
                )
                + sum(C[i, j] * u[j] for j in range(count))
                for i in range(count)
            ]
        boundary = u if boundary is None else boundary
        initial = [v.subs(t, 0) for v in u] if initial is None and u is not None else initial
        if boundary is None:
            raise ValueError("Dirichlet data are required")
        self._functions = {}

        def compile_field(name, values, shape):
            expressions = [expr(v) for v in values]
            self._functions[name] = (
                [sp.lambdify((x, y, z, t), v, "numpy") for v in expressions],
                shape,
            )

        compile_field("a", [a], ())
        compile_field("d", [d], ())
        compile_field("b", b, (dim,))
        compile_field("C", list(C), (count, count))
        compile_field("grad_a", [sp.diff(a, v) for v in coordinates], (dim,))
        for name, values in [("f", source), ("g", boundary), ("initial", initial), ("u", u)]:
            if values is not None:
                if len(values) != count:
                    raise ValueError(f"{name} needs {count} components")
                compile_field(name, values, (count,))
        if u is not None:
            compile_field("gradient", [sp.diff(v, c) for v in u for c in coordinates], (count, dim))

    def _evaluate(self, name, p, t):
        x = p[..., 0]
        y = p[..., 1] if self.dim > 1 else bm.zeros_like(x)
        z = p[..., 2] if self.dim > 2 else bm.zeros_like(x)
        functions, shape = self._functions[name]
        values = [
            bm.broadcast_to(bm.asarray(f(x, y, z, t), dtype=p.dtype), x.shape) for f in functions
        ]
        return bm.stack(values, axis=-1).reshape(x.shape + shape) if shape else values[0]

    @cartesian
    def diffusion_coef(self, p, t=0.0):
        return self._evaluate("a", p, t)

    @cartesian
    def diffusion_gradient(self, p, t=0.0):
        return self._evaluate("grad_a", p, t)

    @cartesian
    def convection_coef(self, p, t=0.0):
        return self._evaluate("b", p, t)

    @cartesian
    def reaction_coef(self, p, t=0.0):
        return self._evaluate("C", p, t)

    @cartesian
    def capacity_coef(self, p, t=0.0):
        return self._evaluate("d", p, t)

    @cartesian
    def source(self, p, t=0.0):
        return self._evaluate("f", p, t)

    @cartesian
    def solution(self, p, t=0.0):
        return self._evaluate("u", p, t)

    @cartesian
    def gradient(self, p, t=0.0):
        return self._evaluate("gradient", p, t)

    @cartesian
    def dirichlet(self, p, t=0.0):
        return self._evaluate("g", p, t)

    @cartesian
    def initial(self, p):
        if "initial" not in self._functions:
            raise ValueError("Transient problems require initial data")
        return self._evaluate("initial", p, 0.0)
