"""A: FEALPy mesh import and orientation-sensitive triangle geometry."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from fealpy.backend import backend_manager as bm
from fealpy.mesh import TriangleMesh

from src.geometry import edge_matrix, signed_simplex_measure, signed_triangle_areas


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    bm.set_backend("numpy")
    mesh = TriangleMesh.from_box([0.0, 1.0, 0.0, 1.0], nx=2, ny=2)
    nodes = bm.to_numpy(mesh.entity("node"))
    cells = bm.to_numpy(mesh.entity("cell"))
    signed = signed_triangle_areas(nodes, cells)
    boundary_count = int(np.count_nonzero(bm.to_numpy(mesh.boundary_node_flag())))

    triangle = np.array([[0.0, 0.0], [2.0, 0.0], [0.0, 1.0]])
    flipped = triangle[[0, 2, 1]]
    degenerate = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])

    result = {
        "backend": "numpy",
        "dtype": str(nodes.dtype),
        "node_count": int(mesh.number_of_nodes()),
        "cell_count": int(mesh.number_of_cells()),
        "boundary_node_count": boundary_count,
        "min_signed_area": float(signed.min()),
        "max_signed_area": float(signed.max()),
        "total_signed_area": float(signed.sum()),
        "manual_edge_matrix": edge_matrix(triangle).tolist(),
        "manual_jacobian_determinant": float(np.linalg.det(edge_matrix(triangle))),
        "manual_signed_area": signed_simplex_measure(triangle),
        "flipped_signed_area": signed_simplex_measure(flipped),
        "degenerate_signed_area": signed_simplex_measure(degenerate),
    }

    assert result["dtype"] == "float64"
    assert np.isclose(result["total_signed_area"], 1.0)
    assert result["min_signed_area"] > 0.0
    assert np.isclose(result["manual_signed_area"], 1.0)
    assert np.isclose(result["flipped_signed_area"], -1.0)
    assert np.isclose(result["degenerate_signed_area"], 0.0)

    output = ROOT / "results" / "a_geometry.json"
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

