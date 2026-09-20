"""Bounded integration checks for the vector solver and its public data API."""

import json, sys, tempfile
from pathlib import Path
import numpy as np
import fealpy, fealpy.fem
from fealpy.backend import backend_manager as bm


def run(output):
    exports = dict(vars(fealpy.fem))
    paths = list(sys.path)
    bm.set_backend("pytorch")
    from cdr import CDRSolver, SymbolicCDR, box_mesh
    from cdr.cases import CoupledData

    assert bm.get_current_backend().backend_name == "pytorch"
    assert sys.path == paths and set(exports) == set(vars(fealpy.fem))
    assert all(vars(fealpy.fem)[k] is v for k, v in exports.items())
    bm.set_backend("numpy")
    records = []
    with tempfile.TemporaryDirectory() as tmp:
        for dim in (1, 2, 3):
            for p in (1, 2):
                data = CoupledData(dim, "quadratic", "linear")
                mesh = box_mesh(data, 3)
                nodes = mesh.node.copy()
                cells = mesh.cell.copy()
                marker = np.arange(mesh.number_of_nodes(), dtype=float)
                mesh.nodedata["u"] = marker
                s = CDRSolver(data, mesh, p=p)
                first = s.solve_time(0.1, 2)[:].copy()
                other = CDRSolver(data, mesh, p=p)
                other.solve_time(0.2, 2)
                np.testing.assert_allclose(s.solve_time(0.1, 2)[:], first, atol=1e-12)
                s.export_vtu(s.uh, Path(tmp) / f"{dim}_{p}.vtu")
                assert mesh.nodedata["u"] is marker
                np.testing.assert_array_equal(mesh.node, nodes)
                np.testing.assert_array_equal(mesh.cell, cells)
                if p == 2:
                    assert s.errors(0.1)[0] < 1e-10
                records.append(dict(dim=dim, p=p, L2=s.errors(0.1)[0]))
        data = SymbolicCDR(
            1,
            components=4,
            a="1+t+x*x",
            b=["1+t"],
            reaction=2,
            source=[1, 2, 3, 4],
            boundary=[0, 0, 0, 0],
            initial=[0, 0, 0, 0],
        )
        s = CDRSolver(data, box_mesh(data, 8))
        s.solve_time(0.1, 4, output=Path(tmp) / "four")
        import vtk

        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(str(Path(tmp) / "four/solution_0004.vtu"))
        reader.Update()
        assert all(
            reader.GetOutput().GetPointData().GetArray(f"u_{i}") is not None for i in range(4)
        )
        assert np.max(s.uh[:]) > 0
        from fealpy.solver import spsolve

        def modifying(A, F):
            values = spsolve(A, F, solver="scipy")
            F[:] = 999
            return values

        s = CDRSolver(data, box_mesh(data, 8), linear_solver=modifying)
        s.solve_time(0.1, 4)
        assert s.history[-1]["residual"] < 1e-10
        bm.set_backend("pytorch")
        try:
            s.solve_time(0.1, 4)
        except RuntimeError:
            pass
        else:
            raise AssertionError("Mixed backend reuse must be rejected")
        assert bm.get_current_backend().backend_name == "pytorch"
        bm.set_backend("numpy")
        # Strong-advection, time-dependent manufactured data also exercises SUPG mass/load.
        for dim in (1, 2, 3):
            names = ["x", "y", "z"][:dim]
            w = "(1+" + "+".join(v + "**2" for v in names) + ")"
            data = SymbolicCDR(
                dim,
                exact=[f"(1+t)*{w}", f"2*(1+t)*{w}"],
                a=".001*(1+t)",
                b=["10*(1+t)"] * dim,
                reaction=[[2, -1], [-1, 2]],
                d="1+t",
            )
            s = CDRSolver(data, box_mesh(data, 3), p=2)
            s.solve_time(0.1, 2)
            assert s.errors(0.1)[0] < 1e-10
        # A time-dependent capacity multiplies u_t, not the derivative of d*u.
        data = SymbolicCDR(1, exact=["exp(t)*x*(1-x)"], d="1+t")
        s = CDRSolver(data, box_mesh(data, 8), p=2)
        s.solve_time(0.2, 40)
        assert s.errors(0.2)[0] < 1e-5
        # Non-unit boxes exercise geometric boundary detection away from 0 and 1.
        for dim in (1, 2, 3):
            names = ["x", "y", "z"][:dim]
            expression = "(1+t)*(1+" + "+".join(v + "**2" for v in names) + ")"
            data = SymbolicCDR(dim, exact=[expression], box=[-1.0, 2.0] * dim)
            s = CDRSolver(data, box_mesh(data, 3), p=2)
            s.solve_time(0.1, 2)
            assert s.errors(0.1)[0] < 1e-10
    result = dict(
        fealpy_path=fealpy.__file__,
        records=records,
        import_isolation=True,
        mesh_isolation=True,
        repeated_interleaved=True,
        callback_isolation=True,
        backend_guard=True,
        four_component_export=True,
        no_exact_solution=True,
        strong_transient_consistency=True,
        capacity_semantics=True,
        nonunit_boxes=True,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    run(Path(sys.argv[1]) if len(sys.argv) > 1 else Path("results_final/interfaces_installed.json"))
