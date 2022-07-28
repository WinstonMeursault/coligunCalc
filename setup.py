from distutils.core import Extension, setup

from Cython.Build import cythonize

setup(
    name = "coilgunCalc",
    ext_modules=cythonize(["basics.pyx", "core.pyx"], annotate=False)
)