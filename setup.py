from distutils.core import Extension, setup

from Cython.Build import cythonize

setup(
    name = "Consolidated Loop",
    ext_modules=cythonize(["basics.pyx", "core.pyx"])
)