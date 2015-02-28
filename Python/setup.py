import sys
sys.path.insert(0, "..")
print sys.path

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

ext_modules = cythonize("**/*.pyx", exclude="sim_*.pyx")

from numpy.distutils.misc_util import get_numpy_include_dirs
numpy_demo = [Extension("*",["sim_free.pyx"],
                        include_dirs=get_numpy_include_dirs())]
ext_modules.extend(cythonize(numpy_demo))

setup(
    name = "My hello app",
    ext_modules = ext_modules, # accepts a glob pattern
)
