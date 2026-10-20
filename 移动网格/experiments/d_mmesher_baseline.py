"""D: smallest practical two-dimensional MMesher/GFMMPDE baseline."""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from fealpy.backend import backend_manager as bm
from fealpy.decorator import cartesian
from fealpy.functionspace import LagrangeFESpace
from fealpy.mesh import TriangleMesh
from fealpy.mmesh import Config, MMesher

from src.geometry import signed_triangle_areas


ROOT = Path(__file__).resolve().parents[1]


@cartesian
def layer(p):
    x, y = p[..., 0], p[..., 1]
    return bm.tanh(40.0 * (1.0 - x - y))


def main() -> None:
    bm.set_backend("numpy")
    mesh = TriangleMesh.from_box_cross_mesh([0.0, 1.0, 0.0, 1.0], nx=4, ny=4)
    mesh.meshdata["vertices"] = bm.array(
        [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]], dtype=bm.float64
    )
    cells = bm.to_numpy(mesh.entity("cell"))
    initial_nodes = bm.to_numpy(mesh.entity("node")).copy()
    initial_areas = signed_triangle_areas(initial_nodes, cells)

    space = LagrangeFESpace(mesh, p=1)
    uh = space.interpolate(layer)
    initial_error = float(mesh.error(layer, uh, q=5))

    config = Config(active_method="GFMMPDE", is_pre=False, maxit=3, t_max=0.1, tau=0.01)
    mover = MMesher(mesh=mesh, uh=uh, space=space, beta=0.5, config=config)
    started = time.perf_counter()
    mover.initialize()
    moved_mesh, moved_uh = mover.instance.mesh_redistributor(maxit=3)
    elapsed = time.perf_counter() - started

    moved_nodes = bm.to_numpy(moved_mesh.entity("node"))
    moved_areas = signed_triangle_areas(moved_nodes, cells)
    moved_space = LagrangeFESpace(moved_mesh, p=1)
    moved_function = moved_space.function()
    moved_function[:] = moved_uh
    moved_error = float(moved_mesh.error(layer, moved_function, q=5))
    result = {
        "method": "GFMMPDE",
        "mesh_nx": 4,
        "mesh_ny": 4,
        "iterations_requested": 3,
        "node_count": int(moved_mesh.number_of_nodes()),
        "cell_count": int(moved_mesh.number_of_cells()),
        "initial_min_signed_area": float(initial_areas.min()),
        "moved_min_signed_area": float(moved_areas.min()),
        "initial_total_signed_area": float(initial_areas.sum()),
        "moved_total_signed_area": float(moved_areas.sum()),
        "maximum_node_displacement": float(np.linalg.norm(moved_nodes - initial_nodes, axis=1).max()),
        "initial_interpolation_l2_error": initial_error,
        "moved_interpolation_l2_error": moved_error,
        "elapsed_seconds": elapsed,
        "discrete_energy_available": False,
        "energy_note": "GFMMPDE baseline does not expose a comparable discrete mesh-energy history.",
    }
    assert result["initial_min_signed_area"] > 0.0
    assert result["moved_min_signed_area"] > 0.0
    assert np.isclose(result["moved_total_signed_area"], 1.0, atol=1e-10)

    (ROOT / "results" / "d_mmesher_baseline.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
