"""Atomic text publication for standalone WGS batch entry points."""
from contextlib import contextmanager
import os
from pathlib import Path
import tempfile


@contextmanager
def atomic_text(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(prefix=path.name+'.', suffix='.tmp', dir=path.parent)
    try:
        with os.fdopen(descriptor, 'w') as handle:
            yield handle
        os.replace(name, path)
    finally:
        Path(name).unlink(missing_ok=True)
