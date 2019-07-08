Fast spectrogram phase recovery using Local Weighted Sums (LWS)
===============================================================

Author: Jonathan Le Roux -- 2008-2019

[![PyPI version](https://badge.fury.io/py/lws.svg)](https://badge.fury.io/py/lws)

The LWS software includes the following files:

    readme.md                    # this file
    LICENSE.txt                  # Apache 2 license file
    lwslib/                      # C/C++ library
      -lwslib.cpp                  # core functions
      -lwslib.h                    # header file
    matlab/                      # Matlab/Mex wrapper
      -run_lws.m                   # Matlab example script
      -build_asymmetric_windows.m  # code to build assymetric windows as in RTISI-LA (Matlab) 
      -create_weights.m            # code to create complex weights used in LWS (Matlab)
      -istft.m                     # inverse STFT code (matlab)
      -stft.m                      # STFT code (matlab)
      -batch_lws.cpp               # mex file for LWS
      -nofuture_lws.cpp            # mex file for "no future" LWS initialization
      -online_lws.cpp              # mex file for online LWS
    python/                      # Python module
      -LICENSE.txt                 # Apache 2 license file
      -MANIFEST.in                 # Manifest file specifying files to be distributed with the Python module
      -Makefile                    # Makefile to manage the Python module compilation and building process
      -README.md                   # Readme file for the Python module
      -_version.py                 # Version file
      -lws.pyx                     # Cython source file
      -lwslib.pxd                  # Cython header file
      -setup.py                    # Main module distribution file

License
-------

Copyright (C) 2008-2019 Jonathan Le Roux
Apache 2.0  (http://www.apache.org/licenses/LICENSE-2.0)

Citing this code
----------------

If you use this code, please cite the following papers.

### LWS ###

Jonathan Le Roux, Hirokazu Kameoka, Nobutaka Ono, Shigeki Sagayama,  
"Fast Signal Reconstruction from Magnitude STFT Spectrogram Based on Spectrogram Consistency,"  
in Proc. International Conference on Digital Audio Effects (DAFx), pp. 397--403, Sep. 2010.

    @InProceedings{LeRoux2010DAFx09,
      author =	 {Jonathan {Le Roux} and Hirokazu Kameoka and Nobutaka Ono and Shigeki Sagayama},
      title =	 {Fast Signal Reconstruction from Magnitude {STFT} Spectrogram Based on Spectrogram Consistency},
      booktitle =	 {Proc. International Conference on Digital Audio Effects (DAFx)},
      year =	 2010,
      pages =	 {397--403},
      month =	 sep
    }


### Online LWS, "No future" LWS ###

Jonathan Le Roux, Hirokazu Kameoka, Nobutaka Ono, Shigeki Sagayama,  
"Phase initialization schemes for faster spectrogram-consistency-based signal reconstruction,"  
in Proc. of ASJ Autumn Meeting, 3-10-3, Sep. 2010.

    @InProceedings{LeRoux2010ASJ09,
      author =	 {Jonathan {Le Roux} and Hirokazu Kameoka and Nobutaka Ono and Shigeki Sagayama},
      title =	 {Phase Initialization Schemes for Faster Spectrogram-Consistency-Based Signal Reconstruction},
      year =	 2010,
      booktitle =	 {Proc. Acoustical Society of Japan Autumn Meeting (ASJ)},
      number =	 {3-10-3},
      month =	 mar
    }

Remark
------

The .cpp files are actually C code with some C99 style comments, but the .cpp extension is needed on Windows for mex to acknowledge the c99 flag (with .c, it is discarded, and -ansi used instead, leading to compilation errors)

Acknowledgements
----------------

The recipe to wrap the LWS C code as a python module was largely inspired by Martin Sosic's post: http://martinsosic.com/development/2016/02/08/wrapping-c-library-as-python-module.html

Installation:
-------------

### Matlab 

1) Compiling mex files

    ```sh
    cd matlab/
    mex -I"../lwslib/" -O CFLAGS="\$CFLAGS -std=c99" -output batch_lws ../lwslib/lwslib.cpp batch_lws.cpp
    mex -I"../lwslib/" -O CFLAGS="\$CFLAGS -std=c99" -output online_lws ../lwslib/lwslib.cpp online_lws.cpp
    mex -I"../lwslib/" -O CFLAGS="\$CFLAGS -std=c99" -output nofuture_lws ../lwslib/lwslib.cpp nofuture_lws.cpp
    ```

2) Usage

Please follow/modify `run_lws.m`.
Three steps are implemented, and they can be turned on/off independently:
  * "no future" LWS: phase initialization using LWS updates that only involve past frames
  * online LWS: phase estimation using online LWS updates, corresponding to a fast time-frequency domain version of RTISI-LA
  * batch LWS: phase estimation using batch LWS updates on the whole spectrogram

### Python

1) The easiest way to install `lws` is via `pip`:

    ```sh
    pip install lws
    ```

2) To compile from source using cython (required if one modifies the code):

    ```sh
    cd python
    LWS_USE_CYTHON=1 make
    ```

3) To compile from source using the pre-generated c source file (which was obtained with cython):

    ```sh
    cd python
    make
    ```
    
4) Alternatively, one can first use cython to create a tarball, which can then be installed with pip:

    ```sh
    cd python
    make sdist
    pip install dist/lws-1.2.3.tar.gz
    ```

**Note:** On Windows, the Microsoft Visual C++ Compiler for your version of Python needs to be installed. See [this page](https://wiki.python.org/moin/WindowsCompilers) for more details.

For usage, please refer to `python/README.md`.