#!/usr/bin/env python

import os
from distutils.core import setup, Extension
from numpy.distutils.misc_util import get_numpy_include_dirs

include_dirs = ['/usr/local/include'] + get_numpy_include_dirs()
swig_opts = ["-c++"] + ["-I" + d for d in include_dirs]

ocrorast = Extension('_ocrorast',
                       library_dirs = ['.'],
                       libraries = ['ocrorast','iulib','m'],
                       swig_opts = swig_opts,
                       include_dirs = include_dirs,
                       sources=['ocrorast.i'])

setup(name = 'ocrorast',
      version = '0.1',
      author      = "Thomas Breuel",
      description = """OCRopus layout analysis""",
      ext_modules = [ocrorast],
      py_modules = ["ocrorast","ocropreproc"],
      )
