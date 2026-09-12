"""Symbolic data for scalar diffusion, vector convection and reaction."""
import numpy as np
import sympy as sp
from fealpy.decorator import cartesian


class CDRData:
    """Data on a box for d*u_t - div(a*grad(u)) + b.grad(u) + c*u = f.

    Expressions use coordinates x0, x1, ... and time t. Supply source and
    boundary data for a physical problem; omit source only with an exact solution.
    Numerical evaluation uses NumPy on CPU.
    """
    def __init__(self, dim, a=1, b=None, c=0, d=1, *, exact=None,
                 source=None, boundary=None, initial=None, box=None):
        self.dim = dim
        self.x = sp.symbols(f'x0:{dim}')
        self.t = sp.Symbol('t')
        self.box = list(box if box is not None else [0, 1]*dim)
        symbols = {str(v): v for v in (*self.x, self.t)}
        def expression(value):
            expr = sp.sympify(value)
            return expr.xreplace({v: symbols[str(v)] for v in expr.free_symbols if str(v) in symbols})
        a, c, d = map(expression, (a, c, d))
        b = [expression(v) for v in (b if b is not None else [0]*dim)]
        if len(b) != dim:
            raise ValueError('Convection needs one component per spatial dimension')
        u = expression(exact) if exact is not None else None
        if source is None:
            if u is None:
                raise ValueError('Supply source or an exact manufactured solution')
            source = d*sp.diff(u, self.t) + c*u + sum(
                -sp.diff(a*sp.diff(u, x), x) + bi*sp.diff(u, x)
                for x, bi in zip(self.x, b))
        if boundary is None:
            boundary = u
        if boundary is None:
            raise ValueError('Dirichlet data are required')
        if initial is None and u is not None:
            initial = u.subs(self.t, 0)
        source, boundary = expression(source), expression(boundary)
        initial = expression(initial) if initial is not None else None
        self.expressions = dict(a=a, b=b, c=c, d=d, source=source, exact=u)
        self.diffusion = self._scalar(a)
        self.reaction = self._scalar(c)
        self.capacity = self._scalar(d)
        self.source = self._scalar(source)
        self.dirichlet = self._scalar(boundary)
        self.convection = self._vector(b)
        self.initial = self._scalar(initial) if initial is not None else None
        self.solution = self._scalar(u) if u is not None else None
        self.gradient = self._vector([sp.diff(u, x) for x in self.x]) if u is not None else None

    def domain(self):
        return list(self.box)

    def _scalar(self, expr):
        fun = sp.lambdify((*self.x, self.t), expr, 'numpy')
        @cartesian
        def evaluate(points, t=0):
            return np.broadcast_to(np.asarray(fun(*(points[..., i] for i in range(self.dim)), t),
                                              dtype=float), points.shape[:-1])
        return evaluate

    def _vector(self, expressions):
        components = [self._scalar(e) for e in expressions]
        @cartesian
        def evaluate(points, t=0):
            return np.stack([f(points, t) for f in components], axis=-1)
        return evaluate

    def at_time(self, function, t):
        @cartesian
        def evaluate(points):
            return function(points, t)
        return evaluate
