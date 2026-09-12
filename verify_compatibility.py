"""Check FEALPy integration and shared-object isolation."""
import json
import sys
import tempfile
from pathlib import Path
import numpy as np
import sympy as sp
import fealpy
import fealpy.fem
from fealpy.backend import backend_manager as bm
from fealpy.functionspace import LagrangeFESpace

exports = dict(vars(fealpy.fem))
path = list(sys.path)
bm.set_backend('pytorch')
from pde import CDRData
from cdr_lfem_solver import CDRLFEMSolver, box_mesh, direct_solve
assert bm.get_current_backend().backend_name == 'pytorch'
assert sys.path == path
assert all(vars(fealpy.fem).get(k) is v for k, v in exports.items())
assert set(vars(fealpy.fem)) == set(exports)
bm.set_backend('numpy')
records = []
with tempfile.TemporaryDirectory() as tmp:
    for dim in (1, 2, 3):
        x = sp.symbols(f'x0:{dim}', real=True)
        t = sp.Symbol('t', real=True)
        u = (1+t)*(1+sum(v*v for v in x))
        data = CDRData(dim, exact=u, d=1+t+sum(x))
        for p in (1, 2):
            mesh = box_mesh(data, 3)
            space = LagrangeFESpace(mesh, p=p)
            original = np.arange(mesh.number_of_nodes(), dtype=float)
            mesh.nodedata['u'] = original
            mesh.nodedata['other'] = original+2
            nodes, cells = mesh.node.copy(), mesh.cell.copy()
            s = CDRLFEMSolver(space, data)
            uh = s.solve_time(.1, 2)
            expected = uh[:].copy()
            s.export_vtu(uh, Path(tmp)/f'{dim}_{p}.vtu')
            assert mesh.nodedata['u'] is original
            np.testing.assert_array_equal(mesh.nodedata['other'], original+2)
            np.testing.assert_array_equal(mesh.node, nodes)
            np.testing.assert_array_equal(mesh.cell, cells)
            other = CDRLFEMSolver(space, CDRData(dim, source=1, boundary=0))
            other.solve()
            np.testing.assert_allclose(s.solve_time(.1, 2)[:], expected, atol=1e-12)
            err = float(mesh.error(data.at_time(data.solution, .1), uh, q=6))
            if p == 2:
                assert err < 1e-11, err
            import vtk
            from vtk.util.numpy_support import vtk_to_numpy
            reader = vtk.vtkXMLUnstructuredGridReader()
            reader.SetFileName(str(Path(tmp)/f'{dim}_{p}.vtu')); reader.Update()
            values = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray('u'))
            target = np.empty(mesh.number_of_nodes())
            target[mesh.cell] = uh(np.eye(dim+1))
            np.testing.assert_allclose(values, target)
            records.append(dict(dim=dim, p=p, L2=err, max_residual=max(v['residual'] for v in s.history)))
    def modifying_solver(A, F):
        result = direct_solve(A, F)
        F[:] = 999
        return result
    s = CDRLFEMSolver(space, data, linear_solver=modifying_solver)
    s.solve_time(.1, 2)
    assert s.history[-1]['residual'] < 1e-10
    bm.set_backend('pytorch')
    try:
        s.solve()
    except RuntimeError as exc:
        assert 'NumPy' in str(exc)
    else:
        raise AssertionError('Backend change must be explicit')
    assert bm.get_current_backend().backend_name == 'pytorch'
    bm.set_backend('numpy')
result = dict(fealpy_path=fealpy.__file__, imports_preserve_backend=True,
              imports_preserve_sys_path=True, fem_exports_unchanged=True,
              mesh_and_existing_data_unchanged=True, repeated_interleaved_solves=True,
              modifying_callback_isolated=True, backend_change_rejected=True,
              symbolic_assumptions_supported=True, vtk_values_checked=True, records=records)
output = Path(sys.argv[1])
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(result, indent=2), encoding='utf-8')
print(json.dumps(result, indent=2))
