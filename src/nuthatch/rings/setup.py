from setuptools import setup
from Cython.Build import cythonize

setup(
    name='cythonpoly',
    ext_modules=cythonize("C:/Users/laz01/nuthatch/src/nuthatch/rings/polynomials.py"),
)