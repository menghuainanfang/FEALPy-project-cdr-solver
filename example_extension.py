"""Solve a transient two-dimensional CDR problem without an exact solution."""
from fealpy.backend import backend_manager as bm
from fealpy.functionspace import LagrangeFESpace
from pde import CDRData
from cdr_lfem_solver import CDRLFEMSolver, box_mesh

def main():
    bm.set_backend('numpy')
    pde = CDRData(2, a='1+x0**2+x1**2', b=['1+x0', '1+x1'], c=2,
                  d='1+x0+x1', source=1, boundary=0, initial=0)
    space = LagrangeFESpace(box_mesh(pde, 8), p=1)
    solver = CDRLFEMSolver(space, pde)
    uh = solver.solve_time(end_time=.2, steps=4)
    print('Final range:', uh[:].min(), uh[:].max())
    print('Final solve:', solver.history[-1])


if __name__ == '__main__':
    main()
