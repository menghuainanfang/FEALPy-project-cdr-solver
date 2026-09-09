"""Integration isolation checks for the library candidate, CPU backends only."""
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np
from fealpy.backend import backend_manager as bm

# Import while PyTorch is selected: the module must not reset it to NumPy.
bm.set_backend('pytorch')
import cdr_lfem_solver_1d as module
import cdr_solver
assert bm.get_current_backend().backend_name == 'pytorch'
from fealpy.mesh import IntervalMesh
from fealpy.functionspace import LagrangeFESpace
from fealpy.solver import spsolve
from fealpy.decorator import cartesian
import fealpy, fealpy.fem

exports_before = set(vars(fealpy.fem))
name='fealpy.fem.cdr_lfem_solver_1d'
assert not (Path(fealpy.fem.__file__).parent/'cdr_lfem_solver_1d.py').exists()
spec=importlib.util.spec_from_file_location(name, Path(__file__).with_name('cdr_lfem_solver_1d.py'))
staged=importlib.util.module_from_spec(spec)
spec.loader.exec_module(staged)
assert set(vars(fealpy.fem)) == exports_before
Solver=staged.CdrLFEMSolver1D
records=[]
for backend in ('numpy','pytorch'):
    bm.set_backend(backend)
    @cartesian
    def exact(points):
        xx=points[...,0]
        return 1+xx+bm.sin(bm.pi*xx)
    @cartesian
    def grad(points):
        return (1+bm.pi*bm.cos(bm.pi*points[...,0]))[...,None]
    @cartesian
    def a(points):return 1+points[...,0]**2
    @cartesian
    def b(points):return (1+points[...,0])[...,None]
    @cartesian
    def c(points):return 2+points[...,0]
    @cartesian
    def f(points):
        xx=points[...,0]
        return ((1+xx**2)*bm.pi**2*bm.sin(bm.pi*xx)
                +(1-xx)*(1+bm.pi*bm.cos(bm.pi*xx))+(2+xx)*exact(points))
    for p in (1,2):
        if backend == 'pytorch' and p == 2:
            mesh2=IntervalMesh.from_interval_domain([0,1],nx=8)
            try: Solver(LagrangeFESpace(mesh2,p=2),1,0,0,2,0)
            except NotImplementedError: pass
            else: raise AssertionError('PyTorch P2 limitation must be explicit')
            records.append(dict(backend=backend,p=p,supported=False,reason='IntervalMesh P2 negative-step slicing'))
            continue
        errors=[]
        for n in (8,16,32,64):
            mesh=IntervalMesh.from_interval_domain([0,1],nx=n)
            space=LagrangeFESpace(mesh,p=p)
            nodes=bm.to_numpy(mesh.node).copy()
            cells=bm.to_numpy(mesh.cell).copy()
            model=Solver(space,a,b,c,f,exact)
            uh=model.solve()
            l2=float(mesh.error(exact,uh,q=8))
            h1=float(mesh.error(grad,uh.grad_value,q=8))
            assert model.relative_residual<1e-10 and model.boundary_error<1e-12
            np.testing.assert_array_equal(nodes,bm.to_numpy(mesh.node))
            np.testing.assert_array_equal(cells,bm.to_numpy(mesh.cell))
            baseline=bm.to_numpy(uh[:]).copy()
            matrix=model.A.to_scipy().copy()
            # Residual matches an independently evaluated SciPy matrix product.
            assert np.linalg.norm(matrix @ baseline-bm.to_numpy(model.F))<1e-8
            other=Solver(space,1,0,0,2,0)
            other.solve()
            repeated=model.solve()
            np.testing.assert_allclose(baseline,bm.to_numpy(repeated[:]),atol=1e-12)
            assert bm.get_current_backend().backend_name==backend
            errors.append((l2,h1))
        orders=np.log2(np.array(errors[-2])/np.array(errors[-1]))
        assert orders[0]>p+.8 and orders[1]>p-.2
        records.append(dict(backend=backend,p=p,L2_order=float(orders[0]),H1_order=float(orders[1]),max_last_error=errors[-1][0]))
    # Injection uses FEALPy sparse A and backend F, isolated from model-owned storage.
    calls=[]
    def callback(A,F):
        calls.append(type(A).__name__)
        solution=spsolve(A,F,solver='scipy')
        F[:]=999
        return solution
    injected=Solver(space,1,0,0,2,0,linear_solver=callback)
    injected.solve()
    assert calls and float(bm.max(injected.F))<999
    # Reject backend switches rather than silently altering the global manager.
    bm.set_backend('numpy' if backend=='pytorch' else 'pytorch')
    try: injected.solve()
    except RuntimeError:pass
    else:raise AssertionError('Backend switch was not rejected')
    bm.set_backend(backend)

result=dict(fealpy_path=fealpy.__file__,import_preserves_backend=True,
    staged_module_import=True,existing_fem_exports_unchanged=True,
    mesh_unchanged=True,repeated_and_interleaved_solve=True,solver_injection_isolated=True,
    backend_switch_rejected=True,convergence=records,
    exclusions=['PyTorch P2','GPU','jax','cupy','model-manager registration','official upstream test suite'])
output=Path(__file__).with_name('results_local' if 'site-packages' not in fealpy.__file__ else 'results')
output.joinpath('integration_verification.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
print('PASS: CPU backend compatibility and integration isolation checks')
