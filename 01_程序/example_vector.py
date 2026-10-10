"""Run a coupled time-dependent CDR example and export a ParaView series."""

import argparse
from pathlib import Path
from fealpy.backend import backend_manager as bm
from cdr import CDRSolver, box_mesh
from cdr.cases import CoupledData


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dim", type=int, choices=(1, 2, 3), default=2)
    parser.add_argument("--n", type=int, default=8)
    parser.add_argument("--p", type=int, choices=(1, 2), default=2)
    parser.add_argument("--steps", type=int, default=20)
    parser.add_argument("--scheme", choices=("bdf2", "backward_euler"), default="bdf2")
    parser.add_argument("--output", type=Path, default=Path("output/vector"))
    args = parser.parse_args()
    bm.set_backend("numpy")
    pde = CoupledData(args.dim)
    solver = CDRSolver(pde, box_mesh(pde, args.n), p=args.p)
    solver.solve_time(0.5, args.steps, scheme=args.scheme, output=args.output)
    print("L2 / H1 seminorm:", solver.errors(0.5))
    print("Last step:", solver.history[-1])
    print("ParaView:", args.output / "solution.pvd")


if __name__ == "__main__":
    main()
