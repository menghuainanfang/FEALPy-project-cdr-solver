"""Vector convection-diffusion-reaction finite elements."""

from .solver import CDRSolver, box_mesh
from .pde import CDRModel, SymbolicCDR

__all__ = ["CDRSolver", "box_mesh", "CDRModel", "SymbolicCDR"]
