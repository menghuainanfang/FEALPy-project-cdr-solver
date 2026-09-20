"""Explicit coordinate-based FEALPy data for a coupled manufactured system."""

from fealpy.backend import backend_manager as bm
from fealpy.decorator import cartesian
from .pde import CDRModel


class LayerData(CDRModel):
    """Stable evaluation of -eps*u_xx+u_x=0 on the unit box."""

    def __init__(self, dim=1, eps=0.001):
        super().__init__(dim, 1)
        self.eps = eps

    @cartesian
    def solution(self, p, t=0.0):
        x = p[..., 0]
        tail = bm.exp(bm.asarray(-1 / self.eps))
        return ((bm.exp((x - 1) / self.eps) - tail) / (1 - tail))[..., None]

    @cartesian
    def gradient(self, p, t=0.0):
        value = bm.zeros(p.shape[:-1] + (1, self.dim), dtype=p.dtype)
        value[..., 0, 0] = bm.exp((p[..., 0] - 1) / self.eps) / (
            self.eps * (-bm.expm1(bm.asarray(-1 / self.eps)))
        )
        return value

    @cartesian
    def diffusion_coef(self, p, t=0.0):
        return bm.full_like(p[..., 0], self.eps)

    @cartesian
    def diffusion_gradient(self, p, t=0.0):
        return bm.zeros_like(p)

    @cartesian
    def convection_coef(self, p, t=0.0):
        value = bm.zeros_like(p)
        value[..., 0] = 1.0
        return value

    @cartesian
    def reaction_coef(self, p, t=0.0):
        return bm.zeros(p.shape[:-1] + (1, 1), dtype=p.dtype)

    @cartesian
    def capacity_coef(self, p, t=0.0):
        return bm.ones_like(p[..., 0])

    @cartesian
    def source(self, p, t=0.0):
        return bm.zeros(p.shape[:-1] + (1,), dtype=p.dtype)


class CoupledData(CDRModel):
    """Two components with time-dependent a,b,C,d in dimensions 1/2/3."""

    def __init__(self, dim=2, profile="sine", temporal="exp"):
        super().__init__(dim, 2)
        self.profile, self.temporal = profile, temporal

    def _spatial(self, p):
        x = p[..., 0]
        coordinates = [x]
        if self.dim > 1:
            y = p[..., 1]
            coordinates.append(y)
        if self.dim > 2:
            z = p[..., 2]
            coordinates.append(z)
        if self.profile == "quadratic":
            w = 1 + sum(v * v for v in coordinates)
            grad = bm.stack([2 * v for v in coordinates], axis=-1)
            lap = bm.full_like(x, 2 * self.dim)
        else:
            product = bm.ones_like(x)
            for v in coordinates:
                product *= bm.sin(bm.pi * v)
            w = 1 + sum(coordinates) + product
            derivatives = []
            for i, v in enumerate(coordinates):
                term = bm.pi * bm.cos(bm.pi * v)
                for j, other in enumerate(coordinates):
                    if j != i:
                        term *= bm.sin(bm.pi * other)
                derivatives.append(1 + term)
            grad = bm.stack(derivatives, axis=-1)
            lap = -self.dim * bm.pi**2 * product
        return w, grad, lap

    def _time(self, t):
        if self.temporal == "steady":
            return 1.0, 0.0
        return (
            (1 + t, 1.0)
            if self.temporal == "linear"
            else (bm.exp(bm.asarray(t)), bm.exp(bm.asarray(t)))
        )

    @cartesian
    def diffusion_coef(self, p, t=0.0):
        return (1 + t) * (1 + bm.sum(p * p, axis=-1))

    @cartesian
    def diffusion_gradient(self, p, t=0.0):
        return 2 * (1 + t) * p

    @cartesian
    def convection_coef(self, p, t=0.0):
        return (1 + t) * (1 + p)

    @cartesian
    def reaction_coef(self, p, t=0.0):
        C = bm.asarray([[2.0, -0.5], [-0.25, 3.0]])
        return bm.broadcast_to((1 + t) * C, p.shape[:-1] + (2, 2))

    @cartesian
    def capacity_coef(self, p, t=0.0):
        return 1 + t + bm.sum(p, axis=-1)

    @cartesian
    def solution(self, p, t=0.0):
        w, _, _ = self._spatial(p)
        s, _ = self._time(t)
        return s * bm.stack([w, 2 * w + 1], axis=-1)

    @cartesian
    def gradient(self, p, t=0.0):
        _, grad, _ = self._spatial(p)
        s, _ = self._time(t)
        return s * bm.stack([grad, 2 * grad], axis=-2)

    @cartesian
    def source(self, p, t=0.0):
        w, grad, lap = self._spatial(p)
        s, ds = self._time(t)
        spatial = -self.diffusion_coef(p, t) * lap
        spatial += bm.sum(
            (self.convection_coef(p, t) - self.diffusion_gradient(p, t)) * grad, axis=-1
        )
        ut = ds * bm.stack([w, 2 * w + 1], axis=-1)
        return (
            self.capacity_coef(p, t)[..., None] * ut
            + s * bm.stack([spatial, 2 * spatial], axis=-1)
            + bm.einsum("...ij,...j->...i", self.reaction_coef(p, t), self.solution(p, t))
        )
