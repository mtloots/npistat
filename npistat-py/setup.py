"""Build the shared npistat C back-end at install time, as an Extension, so pip builds real
platform wheels. The library carries no PyInit and is loaded with ctypes.CDLL by _core.py."""
import glob
import os

from setuptools import Extension, setup

# Paths must be RELATIVE to setup.py: setuptools rejects absolute source paths.
SRC = os.path.join("src", "npistat")

sources = sorted(
    s for s in glob.glob(os.path.join(SRC, "*.c"))
    if os.path.basename(s) != "init.c"
)

# Leading underscore: it carries no PyInit and is loaded with ctypes.CDLL by _core.py.
ext = Extension(
    name="npistat._libnpistat",
    sources=sources,
    include_dirs=[SRC],
    libraries=[] if os.name == "nt" else ["m"],
    optional=False,
)

setup(ext_modules=[ext])
