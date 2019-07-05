from setuptools import setup, Extension
import numpy as np
from codecs import open
import os
import sys

# Fix below for Mac from https://github.com/pandas-dev/pandas/pull/24274
# For mac, ensure extensions are built for macos 10.9 when compiling on a
# 10.9 system or above, overriding distuitls behaviour which is to target
# the version that python was built for. This may be overridden by setting
# MACOSX_DEPLOYMENT_TARGET before calling setup.py
if sys.platform == 'darwin':
    import platform
    from distutils.sysconfig import get_config_var
    from distutils.version import LooseVersion
    if 'MACOSX_DEPLOYMENT_TARGET' not in os.environ:
        current_system = LooseVersion(platform.mac_ver()[0])
        python_target = LooseVersion(
            get_config_var('MACOSX_DEPLOYMENT_TARGET'))
        if python_target < '10.9' and current_system >= '10.9':
            os.environ['MACOSX_DEPLOYMENT_TARGET'] = '10.9'


# Get version number from single source, c.f., https://stackoverflow.com/questions/458550/standard-way-to-embed-version-into-python-package
import re
VERSIONFILE = "_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
        verstr = mo.group(1)
else:
        raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))
            
# This file is largely taken from https://github.com/Martinsos/edlib/blob/a77e81678abd9392f1e13ec8585831721a1f354a/bindings/python/setup.py

# Build directly from cython source file(s) if user wants so (probably for some experiments).
# Otherwise, pre-generated c source file(s) are used.
# User has to set environment variable LWS_USE_CYTHON.
# e.g.: LWS_USE_CYTHON=1 python setup.py install
cmdclass = {}
USE_CYTHON = os.getenv('LWS_USE_CYTHON', False)
if USE_CYTHON:
    from Cython.Distutils import build_ext
    lws_module_src = "lws.pyx"
    cmdclass['build_ext'] = build_ext
else:
    lws_module_src = "lws.bycython.cpp"


long_description = ""
# Load README into long description.
here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    # Information
    name = "lws",
    description = "Fast spectrogram phase reconstruction using Local Weighted Sums",
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    version = verstr,
    url = "https://github.com/Jonathan-LeRoux/lws",
    download_url = "https://github.com/Jonathan-LeRoux/lws/archive/{}.tar.gz".format(verstr),
    author = "Jonathan Le Roux",
    author_email = "leroux@merl.com",
    license = "Apache 2.0",
    keywords = ['phase', 'reconstruction', 'stft', 'short-term Fourier Transform', 'spectrogram'],
    # Build instructions
    ext_modules = [Extension("lws",
                             sources=[lws_module_src,"lwslib/lwslib.cpp"],
                             include_dirs=["lwslib/",np.get_include()],
                             language="c++",
                             extra_compile_args=["-O3"])],
    cmdclass = cmdclass,    
    classifiers=[
        "Programming Language :: Python",
        "Intended Audience :: Developers",
        "Topic :: Multimedia :: Sound/Audio :: Analysis",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ]
)
