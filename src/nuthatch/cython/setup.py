from setuptools import setup
from Cython.Build import cythonize

setup(
    name='cpoly',
    package_dir={'nuthatch/src/nuthatch/cython': ''},
    ext_modules=cythonize("C:/Users/laz01/nuthatch/src/nuthatch/cython/cpoly.pyx"),
)