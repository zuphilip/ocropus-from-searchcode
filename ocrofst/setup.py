import os
from distutils.core import setup, Extension
from numpy.distutils.misc_util import get_numpy_include_dirs

baselibs = ['tiff','png','jpeg','SDL','SDL_gfx','SDL_image','m']

include_dirs = ['/usr/local/include'] + get_numpy_include_dirs()
swig_opts = ["-c++"] + ["-I" + d for d in include_dirs]
swiglib = os.popen("swig -swiglib").read()[:-1]

ocrofstll = Extension('_ocrofstll',
        libraries = ['iulib','ocrofst']+baselibs,
        swig_opts = swig_opts,
        include_dirs = include_dirs,
        extra_link_args = ["-L."],
        sources=['ocrofstll.i'])

setup (name = 'ocrofstll',
       version = '0.4',
       author      = "Thomas Breuel",
       description = """ocrofst library bindings""",
       ext_modules = [ocrofstll],
       scripts = ["ocrofst-search"],
       py_modules = ["ocrofstll","ocrofst"],
       )
