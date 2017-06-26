LWS
===

**Fast spectrogram phase recovery using Local Weighted Sums (LWS)**

Author: Jonathan Le Roux -- 2008-2017

| LWS is a C/C++ library for which this package is a Python wrapper.
| A Matlab/Mex wrapper is also available.

License
-------

| Copyright (C) 2008-2017 Jonathan Le Roux
| Apache 2.0  (http://www.apache.org/licenses/LICENSE-2.0)

Reference
---------

If you use this code, please cite the following papers.

Batch LWS
~~~~~~~~~

| Jonathan Le Roux, Hirokazu Kameoka, Nobutaka Ono, Shigeki Sagayama, 
| "Fast Signal Reconstruction from Magnitude STFT Spectrogram Based on Spectrogram Consistency," 
| in Proc. International Conference on Digital Audio Effects (DAFx), pp. 397--403, Sep. 2010.

.. code::

    @InProceedings{LeRoux2010DAFx09,
      author =	 {Jonathan {Le Roux} and Hirokazu Kameoka and Nobutaka Ono and Shigeki Sagayama},
      title =	 {Fast Signal Reconstruction from Magnitude {STFT} Spectrogram Based on Spectrogram Consistency},
      booktitle =	 {Proc. International Conference on Digital Audio Effects (DAFx)},
      year =	 2010,
      pages =	 {397--403},
      month =	 sep
    }
    


Online LWS, "No future" LWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~

| Jonathan Le Roux, Hirokazu Kameoka, Nobutaka Ono, Shigeki Sagayama, 
| "Phase initialization schemes for faster spectrogram-consistency-based signal reconstruction," 
| in Proc. of ASJ Autumn Meeting, 3-10-3, Sep. 2010.

.. code::

    @InProceedings{LeRoux2010ASJ09,
      author =	 {Jonathan {Le Roux} and Hirokazu Kameoka and Nobutaka Ono and Shigeki Sagayama},
      title =	 {Phase Initialization Schemes for Faster Spectrogram-Consistency-Based Signal Reconstruction},
      year =	 2010,
      booktitle =	 {Proceedings of the Acoustical Society of Japan Autumn Meeting (ASJ)},
      number =	 {3-10-3},
      month =	 mar
    }
    
Installation
------------

1. The easiest way to install `lws` is via `pip`:  

        .. code-block:: bash

            pip install lws

2. To compile from source using cython (required if one modifies the code):  

        .. code-block:: bash

            cd python
            LWS_USE_CYTHON=1 make 

3. To compile from source using the pre-generated c source file (which was obtained with cython): 

        .. code-block:: bash

            cd python
            make

4. Alternatively, one can first use cython to create a tarball, which can then be installed by pip:  

        .. code-block:: bash

            cd python
            make sdist
            pip install dist/lws-1.0.tar.gz


Usage
-----
.. code:: python

    import lws
    import numpy as np
    
    lws_processor=lws.lws(512,128, mode="speech") # 512: window length; 128: window shift
    X = lws_processor.stft(x) # where x is a single-channel waveform
    X0 = np.abs(X) # Magnitude spectrogram
    print('{:6}: {:5.2f} dB'.format('Abs(X)', lws_processor.get_consistency(X0))
    X1 = lws_processor.run_lws(X0) # reconstruction from magnitude (in general, one can reconstruct from an initial complex spectrogram)
    print('{:6}: {:5.2f} dB'.format('LWS', lws_processor.get_consistency(X1)))

Options
-------

.. code:: python

    lws_processor=lws.lws(awin_or_fsize, fshift, L = 5, swin = None, look_ahead = 3,
              nofuture_iterations = 0, nofuture_alpha = 1, nofuture_beta = 0.1, nofuture_gamma = 1,
              online_iterations = 0, online_alpha = 1, online_beta = 0.1, online_gamma = 1,
              batch_iterations = 100, batch_alpha = 100, batch_beta = 0.1, batch_gamma = 1,
              symmetric_win = True, mode= None, stft_opts = {})

* ``awin_or_fsize``: either the analysis window, or a window length (in which case the sqrt(hann) window is used)
* ``fshift``: window shift
* ``L``: approximation order in the phase reconstruction algorithm, 5 should be good.
* ``swin``: synthesis window (if None, it gets computed from the analysis window for perfect reconstruction)
* ``look_ahead``: number of look-ahead frames in RTISI-LA-like algorithm, 3 should be good.
* ``xxx_iterations``, ``xxx_alpha``, ``xxx_beta``, ``xxx_gamma``: number of iterations of algorithm xxx (where xxx is one of ``nofuture``, ``online``, or ``batch``), and parameters alpha/beta/gamma of the decreasing sparsity curve that is used to determine which bins get updated at each iteration. Any bin with magnitude larger than a given threshold is updated, others are ignored (``thresholds = alpha * np.exp(- beta * np.arange(iterations)**gamma)``)
* ``symmetric_win``: determines whether to use a symmetric hann window or not
* ``mode``: ``None``, ``'speech'``, or ``'music'``. This sets default numbers of iterations of each algorithm that seem to be good for speech and music signals. Disclaimer: your mileage may vary.
* ``stft_opts``: ``{'perfectrec':True,'fftsize':self.fsize}``. ``perfectrec``: whether to pad with zeros on each side to ensure perfect reconstruction at the boundaries too. ``fftsize``: can be set longer than frame size to do 0-padding in the FFT.

Three steps are implemented, and they can be turned on/off independently by appropriately setting the corresponding number of iterations:

* "no future" LWS: phase initialization using LWS updates that only involve past frames
* online LWS: phase estimation using online LWS updates, corresponding to a fast time-frequency domain version of RTISI-LA
* LWS: phase estimation using batch LWS updates on the whole spectrogram


Remarks
-------

1) The .cpp files are actually C code with some C99 style comments, but the .cpp extension is needed on Windows for mex to acknowledge the c99 flag (with .c, it is discarded, and -ansi used instead, leading to compilation errors)

2) Because the module is a C extension, it cannot be reloaded (see <http://bugs.python.org/issue1144263>). In Jupyter Notebooks, in particular, autoreload will not work, and the kernel has to be restarted.

Acknowledgements
----------------

The recipe to wrap the LWS C code as a python module was largely inspired by the following post by Martin Sosic: <http://martinsosic.com/development/2016/02/08/wrapping-c-library-as-python-module.html>
