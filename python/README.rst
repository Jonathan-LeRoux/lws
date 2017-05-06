===
LWS
===
Fast spectrogram phase recovery using Local Weighted Sums (LWS)

Author: Jonathan Le Roux -- 2008-2017

LWS is a C/C++ library for which this package is a Python wrapper.

-------
License
-------

Copyright (C) 2008-2017 Jonathan Le Roux
Apache 2.0  (http://www.apache.org/licenses/LICENSE-2.0)

---------
Reference
---------

If you use this code, please cite:
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

------------
Installation
------------
::

    pip install lws

-----
Usage
-----
.. code:: python

    import lws


Please follow/modify run_lws.m.
Three steps are implemented, and they can be turned on/off independently:
  * "no future" LWS: phase initialization using LWS updates that only involve past frames
  * online LWS: phase estimation using online LWS updates, corresponding to a fast time-frequency domain version of RTISI-LA
  * LWS: phase estimation using batch LWS updates on the whole spectrogram


Remark: 

the .cpp files are actually C code with some C99 style comments, but the .cpp extension is needed on Windows for mex to acknowledge the c99 flag (with .c, it is discarded, and -ansi used instead, leading to compilation errors)


