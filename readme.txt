Phase recovery using Local Weighted Sums (LWS)
Author: Jonathan Le Roux -- 2008-2017

The LWS software includes the following files:
build_asymmetric_windows.m  # code to build assymetric windows as in RTISI-LA (Matlab) 
run_lws.m                   # example script (Matlab)
lib/create_weights.m        # code to create complex weights used in LWS (Matlab)
lib/istft.m                 # inverse STFT code (matlab)
lib/lws.c                   # mex file for LWS
lib/lws_functions.c         # core functions
lib/lws_functions.h  
lib/nofuture_lws.c          # mex file for "no future" LWS initialization
lib/online_lws.c            # mex file for online LWS
lib/stft.m                  # STFT code (matlab)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2008-2017 Jonathan Le Roux
%   Apache 2.0  (http://www.apache.org/licenses/LICENSE-2.0) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
---------------------


Installation:

1) Compiling mex files

mex -O CFLAGS="\$CFLAGS -std=c99" -output lws lib/lws_functions.c lib/lws.c
mex -O CFLAGS="\$CFLAGS -std=c99" -output online_lws lib/lws_functions.c lib/online_lws.c
mex -O CFLAGS="\$CFLAGS -std=c99" -output nofuture_lws lib/lws_functions.c lib/nofuture_lws.c

2) Usage

Please follow/modify run_lws.m.
Three steps are implemented, and they can be turned on/off independently:
- "no future" LWS: phase initialization using LWS updates that only involve past frames
- online LWS: phase estimation using online LWS updates, corresponding to a fast time-frequency domain version of RTISI-LA
- LWS: phase estimation using batch LWS updates on the whole spectrogram

