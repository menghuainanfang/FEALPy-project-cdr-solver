"""Reproduce spatial, temporal, coupled and convection-dominated validation."""

import argparse, csv, json, platform, tempfile
from pathlib import Path
from importlib.metadata import version
import numpy as np
import sympy as sp
import fealpy
from fealpy.backend import backend_manager as bm
from cdr import CDRSolver, SymbolicCDR, box_mesh
from cdr.cases import CoupledData, LayerData


def verify_history(solver):
    assert max(v["residual"] for v in solver.history) < 1e-10
    assert max(v["boundary_error"] for v in solver.history) < 1e-11


def run(output):
    output.mkdir(parents=True, exist_ok=True)
    bm.set_backend("numpy")
    spatial, temporal, polynomial, strong = [], [], [], []
    for dim in (1, 2, 3):
        grids = (8, 16, 32) if dim == 1 else ((4, 8, 16) if dim == 2 else (2, 4, 8))
        for stabilization in ("galerkin", "supg"):
            for p in (1, 2):
                previous = None
                for n in grids:
                    data = CoupledData(dim, temporal="steady")
                    solver = CDRSolver(data, box_mesh(data, n), p=p, stabilization=stabilization)
                    solver.solve_steady(t=0.3)
                    l2, h1 = solver.errors(0.3)
                    orders = np.log2(np.array(previous) / [l2, h1]) if previous else [0.0, 0.0]
                    row = dict(
                        dim=dim,
                        p=p,
                        n=n,
                        stabilization=stabilization,
                        dofs=solver.space.number_of_global_dofs(),
                        L2=l2,
                        H1=h1,
                        L2_order=float(orders[0]),
                        H1_order=float(orders[1]),
                        **solver.history[-1],
                    )
                    spatial.append(row)
                    previous = (l2, h1)
                    verify_history(solver)
                    print("space", row, flush=True)
                assert orders[0] > p + 0.6 and orders[1] > p - 0.25, row
        for stabilization in ("galerkin", "supg"):
            data = CoupledData(dim, "quadratic", "linear")
            solver = CDRSolver(data, box_mesh(data, 3), p=2, stabilization=stabilization)
            solver.solve_time(0.2, 4)
            err = solver.errors(0.2)
            assert err[0] < 1e-10 and err[1] < 1e-9, err
            polynomial.append(dict(dim=dim, stabilization=stabilization, L2=err[0], H1=err[1]))
        for scheme in ("backward_euler", "bdf2"):
            previous = None
            for steps in (10, 20, 40, 80):
                data = CoupledData(dim, "quadratic", "exp")
                solver = CDRSolver(data, box_mesh(data, 3), p=2)
                solver.solve_time(0.5, steps, scheme=scheme)
                l2, h1 = solver.errors(0.5)
                order = float(np.log2(previous / l2)) if previous else 0.0
                previous = l2
                verify_history(solver)
                row = dict(
                    dim=dim,
                    scheme=scheme,
                    steps=steps,
                    dt=0.5 / steps,
                    L2=l2,
                    H1=h1,
                    order=order,
                    max_residual=max(v["residual"] for v in solver.history),
                )
                temporal.append(row)
                print("time", row, flush=True)
            expected = 1 if scheme == "backward_euler" else 2
            assert expected - 0.2 < order < expected + 0.2, row
    # Symbolic differentiation independently checks the explicit PDE source.
    points = np.array([[0.13, 0.27, 0.41], [0.32, 0.54, 0.61]])
    x, y, z, t = sp.symbols("x y z t", real=True)
    w = 1 + x * x + y * y + z * z
    symbolic = SymbolicCDR(
        3,
        exact=[sp.exp(t) * w, sp.exp(t) * (2 * w + 1)],
        a=(1 + t) * w,
        b=[(1 + t) * (1 + v) for v in (x, y, z)],
        d=1 + t + x + y + z,
        reaction=[[(1 + t) * 2, -(1 + t) / 2], [-(1 + t) / 4, 3 * (1 + t)]],
    )
    explicit = CoupledData(3, "quadratic")
    np.testing.assert_allclose(
        symbolic.source(points, 0.37), explicit.source(points, 0.37), rtol=1e-13
    )
    # Real off-diagonal coupling: changing C with fixed f changes the result.
    common = dict(exact=["1+x*x", "2+2*x*x"], a=1, b=[1], source=["2+2*x+4*x*x", "7+4*x+11*x*x"])
    coupled = SymbolicCDR(1, reaction=[[2, 1], [1, 5]], **common)
    diagonal = SymbolicCDR(1, reaction=[[2, 0], [0, 5]], **common)
    s1 = CDRSolver(coupled, box_mesh(coupled, 12), p=2)
    s2 = CDRSolver(diagonal, box_mesh(diagonal, 12), p=2)
    u1 = s1.solve_steady()[:].copy()
    u2 = s2.solve_steady()[:].copy()
    assert s1.errors()[0] < 1e-11
    coupling_effect = float(np.max(np.abs(u1 - u2)))
    assert coupling_effect > 1e-2
    # Old layer case plus its full-dimensional extensions; no post-solve clipping.
    for dim in (1, 2, 3):
        eps = 0.001
        data = LayerData(dim, eps)
        n = 32 if dim < 3 else 8
        for stabilization in ("galerkin", "supg"):
            solver = CDRSolver(data, box_mesh(data, n), p=1, stabilization=stabilization)
            uh = solver.solve_steady()
            verify_history(solver)
            row = dict(
                dim=dim,
                n=n,
                stabilization=stabilization,
                minimum=float(uh[:].min()),
                maximum=float(uh[:].max()),
                L2=solver.errors()[0],
                **solver.history[-1],
            )
            strong.append(row)
            print("layer", row, flush=True)
            if dim == 1:
                points = solver.scalar_space.interpolation_points()[:, 0]
                np.savez(output / f"layer_{stabilization}.npz", x=points, u=uh[:])
                if stabilization == "galerkin":
                    assert row["minimum"] < -0.8
                else:
                    assert row["minimum"] >= -1e-10 and row["maximum"] <= 1 + 1e-10
                    exact = data.solution(points[:, None])[:, 0]
                    np.testing.assert_allclose(uh[:], exact, atol=1e-10)
            if stabilization == "supg":
                assert row["minimum"] >= -1e-8 and row["maximum"] <= 1 + 1e-8, row
    # Finite layer thickness: convergence after the layer is progressively resolved.
    resolved = []
    for n in (32, 64, 128):
        data = LayerData(1, 0.02)
        s = CDRSolver(data, box_mesh(data, n), p=1)
        s.solve_steady()
        l2, h1 = s.errors()
        resolved.append(dict(n=n, L2=l2, H1=h1))
    assert resolved[-1]["L2"] < resolved[0]["L2"] / 8
    # Time-dependent coupled SUPG solution for ParaView.
    data = CoupledData(2)
    solver = CDRSolver(data, box_mesh(data, 12), p=2)
    original = np.arange(solver.mesh.number_of_nodes(), dtype=float)
    solver.mesh.nodedata["existing"] = original
    uh = solver.solve_time(0.5, 20, output=output / "vtu")
    verify_history(solver)
    assert solver.mesh.nodedata["existing"] is original and "u" not in solver.mesh.nodedata
    np.savez(
        output / "field.npz",
        node=solver.mesh.node,
        cell=solver.mesh.cell,
        exact=data.solution(solver.mesh.node, 0.5),
        initial=data.initial(solver.mesh.node),
    )
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(output / "vtu/solution_0020.vtu"))
    reader.Update()
    grid = reader.GetOutput()
    actual = vtk_to_numpy(grid.GetPointData().GetArray("u"))
    expected = np.empty((solver.mesh.number_of_nodes(), 2))
    expected[solver.mesh.cell] = uh(np.eye(3))
    np.testing.assert_allclose(actual[:, :2], expected)
    assert actual.shape[1] == 3 and np.all(actual[:, 2] == 0)
    np.save(output / "field_values.npy", expected)
    # The one-component interface also reproduces the scalar Poisson limit.
    scalar = SymbolicCDR(1, exact=["x*(1-x)"], source=[2])
    scalar_solver = CDRSolver(scalar, box_mesh(scalar, 8), p=2)
    scalar_solver.solve_steady()
    assert scalar_solver.errors()[0] < 1e-12
    records = dict(
        spatial=spatial,
        temporal=temporal,
        polynomial=polynomial,
        strong=strong,
        resolved_layer=resolved,
        coupling_effect=coupling_effect,
        source_agreement=True,
        vtu_vector_values=True,
        mesh_data_unchanged=True,
        environment=dict(
            python=platform.python_version(),
            fealpy_path=fealpy.__file__,
            versions={v: version(v) for v in ("fealpy", "numpy", "scipy", "sympy", "vtk")},
        ),
    )
    (output / "verification.json").write_text(json.dumps(records, indent=2), encoding="utf-8")
    for name, rows in [("spatial", spatial), ("temporal", temporal), ("strong_convection", strong)]:
        with (output / f"{name}.csv").open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0])
            writer.writeheader()
            writer.writerows(rows)
    print(
        "PASS: vector space/time orders, SUPG residual consistency, coupling, layers and VTU",
        flush=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(__file__).parent / "results_final")
    run(parser.parse_args().output)
