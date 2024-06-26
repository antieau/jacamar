from setuptools import setup
from Cython.Build import cythonize

setup(
    name='polynomials',
    ext_modules=cythonize("src/nuthatch/rings/polynomials.py"),
)
