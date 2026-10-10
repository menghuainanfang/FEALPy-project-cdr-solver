"""Consistent cell residuals for streamline-upwind Petrov-Galerkin FEM."""

from fealpy.backend import backend_manager as bm


def supg_terms(space, pde, t, q):
    """Return local (spatial blocks, time mass, loads) using tau*b.grad(v)."""
    mesh = space.mesh
    bcs, w = mesh.quadrature_formula(q, "cell").get_quadrature_points_and_weights()
    points = mesh.bc_to_point(bcs)
    measure = mesh.entity_measure("cell")
    grad = space.grad_basis(bcs)
    phi = bm.broadcast_to(space.basis(bcs), grad.shape[:-1])
    a, b, C, d = (
        pde.diffusion_coef(points, t),
        pde.convection_coef(points, t),
        pde.reaction_coef(points, t),
        pde.capacity_coef(points, t),
    )
    speed = bm.linalg.norm(b, axis=-1)
    # Streamline length uses linear barycentric gradients, independent of p.
    linear_grad = mesh.grad_lambda()
    directional = bm.einsum("cqd,cid->cqi", b, linear_grad)
    denominator = bm.sum(bm.abs(directional), axis=-1)
    active = speed > 1e-14
    h = 2 * speed / bm.where(denominator > 0, denominator, 1.0) / space.p
    pe = speed * h / (2 * a)
    small = pe < 1e-3
    safe_pe = bm.where(small, 1.0, pe)
    langevin = bm.where(small, pe / 3 - pe**3 / 45, 1 / bm.tanh(safe_pe) - 1 / safe_pe)
    tau = bm.where(active, h * langevin / (2 * bm.where(active, speed, 1.0)), 0.0)
    streamline = bm.einsum("cqd,cqid->cqi", b, grad)
    test = tau[..., None] * streamline
    laplace = bm.einsum("cqidd->cqi", space.hess_basis(bcs)) if space.p == 2 else bm.zeros_like(phi)
    residual = streamline - a[..., None] * laplace
    residual -= bm.einsum("cqd,cqid->cqi", pde.diffusion_gradient(points, t), grad)

    def integrate(trial):
        return bm.einsum("q,c,cqi,cqj->cij", w, measure, test, trial)

    m = pde.components
    blocks = [
        [integrate(C[..., i, j, None] * phi + (residual if i == j else 0.0)) for j in range(m)]
        for i in range(m)
    ]
    mass = integrate(d[..., None] * phi)
    source = pde.source(points, t)
    loads = [bm.einsum("q,c,cqi,cq->ci", w, measure, test, source[..., i]) for i in range(m)]
    return blocks, mass, loads
