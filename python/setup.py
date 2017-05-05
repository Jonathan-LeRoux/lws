from setuptools import setup, Extension
import numpy as np

setup(
    # Information
    name = "lws",
    version = "1.0.0",
    url = "",
    author = "Jonathan Le Roux",
    license = "Apache 2.0",
    keywords = "phase reconstruction stft",
    # Build instructions
    ext_modules = [Extension("lws",
                             sources=["c/lws_functions.cpp","lws.pyx"],
                             include_dirs=["c/",np.get_include()],
                             language="c++"])
)