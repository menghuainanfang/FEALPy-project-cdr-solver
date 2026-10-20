"""Package the two offline viewers, numerical code and FEALPy source snapshot."""

from datetime import datetime
from importlib.util import find_spec
from pathlib import Path
import shutil
import tempfile
import zipfile
import subprocess


ROOT = Path(__file__).resolve().parents[1]
NAME = "moving_mesh_demos_2026-10-17"


def main():
    fealpy_root = Path(find_spec("fealpy").origin).resolve().parents[1]
    with tempfile.TemporaryDirectory(prefix="moving-mesh-package-") as temporary:
        package = Path(temporary) / NAME
        files = [
            "src/__init__.py", "src/geometry.py", "src/variational_mesh.py",
            "src/poisson_mesh_comparison.py", "src/fealpy_mesh_evolution.py", "experiments/__init__.py",
            "experiments/e_variational_mesh_visual.py", "experiments/mesh_viewer.html",
            "experiments/f_poisson_mesh_comparison.py", "experiments/poisson_viewer.html",
            "tests/test_geometry.py", "tests/test_variational_mesh.py",
            "tests/test_poisson_mesh_comparison.py", "tests/test_fealpy_mesh_evolution.py", "docs/variational_mesh_visual.md",
            "docs/poisson_mesh_comparison.md", "requirements-lock-windows-py313.txt",
        ]
        for prefix, suffixes in [("e_variational_mesh", ["html", "json", "png"]),
                                  ("f_poisson_comparison", ["html", "json", "csv", "png"])]:
            files.extend(f"results/{prefix}.{suffix}" for suffix in suffixes)
        for relative in files:
            target = package / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ROOT / relative, target)
        assets = Path(__file__).with_name("demo_package")
        for source in assets.iterdir():
            shutil.copy2(source, package / source.name)
        vendor = package / "vendor" / "fealpy_source"
        vendor.mkdir(parents=True)
        for directory in ["fealpy", "external_deps"]:
            shutil.copytree(fealpy_root / directory, vendor / directory,
                            ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        for filename in ["setup.py", "pyproject.toml", "requirements.txt", "MANIFEST.in",
                         "README.md", "COPYRIGHT.txt", "GPL.pdf"]:
            shutil.copy2(fealpy_root / filename, vendor / filename)
        revision = subprocess.check_output(["git", "-C", str(fealpy_root), "rev-parse", "HEAD"], text=True).strip()
        (package / "PACKAGE_INFO.txt").write_text(
            f"Packaged: {datetime.now().astimezone().isoformat()}\n"
            f"FEALPy source base: {revision}\n"
            "Includes the local FEALPy source snapshot used by the demonstrations.\n"
            "No Git history, Python environment, paper PDFs or personal study notes included.\n",
            encoding="utf-8")
        archive = ROOT.parent / f"{NAME}.zip"
        with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as output:
            for source in sorted(package.rglob("*")):
                if source.is_file():
                    output.write(source, source.relative_to(Path(temporary)))
        print(f"Archive: {archive}\nSize: {archive.stat().st_size / 1024**2:.2f} MiB")


if __name__ == "__main__":
    main()
