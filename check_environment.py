"""Report interpreter/dependencies and fail early when the project runtime drifts."""
import importlib
import importlib.metadata as metadata
import sys
from pathlib import Path


def check():
    versions = {}
    for line in (Path(__file__).parent / 'requirements.txt').read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        name, expected = line.split('==')
        module = importlib.import_module(name)
        actual = metadata.version(name)
        if actual != expected:
            raise RuntimeError(f'{name}: expected {expected}, found {actual}')
        versions[name] = actual
    print('Python:', sys.executable)
    print('Dependencies:', ', '.join(f'{k} {v}' for k, v in versions.items()))
    print('FEALPy:', importlib.import_module('fealpy').__file__)
    return versions


if __name__ == '__main__':
    check()
