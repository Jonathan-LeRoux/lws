Phase recovery using Local Weighted Sums (LWS)
Author: Jonathan Le Roux -- 2008-2017

The LWS software includes the following files:
readme.txt                   # this file
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
  -README.rst                  # Readme file for the Python module
  -lws.pyx                     # Cython source file
  -lwslib.pxd                  # Cython header file
  -setup.py                    # Main module distribution file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2008-2017 Jonathan Le Roux                   % 
%   Apache 2.0  (http://www.apache.org/licenses/LICENSE-2.0)   % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Citing this code

If you use this code, please cite the following papers.

## LWS ##

Jonathan Le Roux, Hirokazu Kameoka, Nobutaka Ono and Shigeki Sagayama, 
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


## Online LWS, "No future" LWS ##

Jonathan Le Roux, Hirokazu Kameoka, Nobutaka Ono and Shigeki Sagayama, 
"Phase initialization schemes for faster spectrogram-consistency-based signal reconstruction," 
Proc. of ASJ Autumn Meeting, 3-10-3, Sep. 2010.

    @InProceedings{LeRoux2010ASJ09,
      author =	 {Jonathan {Le Roux} and Hirokazu Kameoka and Nobutaka Ono and Shigeki Sagayama},
      title =	 {Phase Initialization Schemes for Faster Spectrogram-Consistency-Based Signal Reconstruction},
      year =	 2010,
      booktitle =	 {Proceedings of the Acoustical Society of Japan Autumn Meeting (ASJ)},
      number =	 {3-10-3},
      month =	 mar
    }

# Remark

The .cpp files are actually C code with some C99 style comments, but the .cpp extension is needed on Windows for mex to acknowledge the c99 flag (with .c, it is discarded, and -ansi used instead, leading to compilation errors)

# Acknowledgements

The recipe to wrap the LWS C code as a python module was largely inspired by Martin Sosic's post: http://martinsosic.com/development/2016/02/08/wrapping-c-library-as-python-module.html

# Installation:

## Matlab 
1) Compiling mex files

    cd matlab/
    mex -I"../lwslib/" -O CFLAGS="\$CFLAGS -std=c99" -output batch_lws ../lwslib/lwslib.cpp batch_lws.cpp
    mex -I"../lwslib/" -O CFLAGS="\$CFLAGS -std=c99" -output online_lws ../lwslib/lwslib.cpp online_lws.cpp
    mex -I"../lwslib/" -O CFLAGS="\$CFLAGS -std=c99" -output nofuture_lws ../lwslib/lwslib.cpp nofuture_lws.cpp

2) Usage

Please follow/modify run_lws.m.
Three steps are implemented, and they can be turned on/off independently:
  * "no future" LWS: phase initialization using LWS updates that only involve past frames
  * online LWS: phase estimation using online LWS updates, corresponding to a fast time-frequency domain version of RTISI-LA
  * batch LWS: phase estimation using batch LWS updates on the whole spectrogram

## Python

1) Compiling using cython:
    cd python
    make

2) Alternatively, one can install from the tarball without needing cython:
    pip install dist/lws-1.0.tar.gz

3) Install using pip
    pip install lws
