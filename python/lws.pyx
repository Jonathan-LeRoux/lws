#cython: language_level=3
from __future__ import division
cimport lwslib
cimport numpy as np
import numpy as np
import scipy

__version__ = "1.2.6"

def hann(n,symmetric=True,use_offset = False):
    if symmetric:
        w = 0.5*(1 - np.cos(2*np.pi*(np.arange(1,2*n,2))/(2*n)))
    else:
        if use_offset:
            offset = 1
        else:
            offset = 0
        w = 0.5*(1 - np.cos(2*np.pi*(np.arange(n)+offset)/n))
    return w


def synthwin(awin,fshift,swin=None):
    # returns the normalized synthesis window for perfect reconstruction
    fsize = len(awin)
    Q = np.int(np.ceil(np.float(fsize) / np.float(fshift)))
    if swin is None:
        swin = awin
    twin = awin * swin                              
    tmp_fsize = Q*fshift
    
    w=np.hstack([twin,np.zeros((tmp_fsize-fsize,))])
    w=np.sum(np.reshape(w,(Q,fshift)),axis=0)
    w=np.tile(w,(1,Q))[0,:fsize]
    
    if min(w) <= 0:
        raise ValueError('The normalizer is not strictly positive')

    swin = swin / w

    return swin


def stft(x,fsize,fshift,awin,fftsize=None,perfectrec=False):
    # STFT with a fixed frame shift

    if len(np.shape(x)) != 1: # no multi-channel input
        raise ValueError('We only deal with single channel signals here')

    if fftsize is None:
        fftsize = fsize
    if fftsize % 2 ==1:
        raise ValueError('Odd ffts not supported.')

    if perfectrec is True:
        # pad beginning so that original signal starts at a frame start
        residual_size = fsize % fshift
        if residual_size == 0:
            pre_pad = np.zeros((fsize - fshift,))
        else:
            pre_pad = np.zeros((fsize - residual_size,))
        # pad end just enough to complete the last frame including original signal;
        # the next frame, if added, would not contain any signal
        if len(x) % fshift == 0:
            post_pad = np.zeros((0,))
        else:
            post_pad = np.zeros((fshift - len(x) % fshift,))
        x = np.hstack((pre_pad, x, post_pad))
        # after padding, len(x) is a multiple of fshift
        M = len(x)//fshift
    else:
        # Do not care about all signal being in regions where frames give complete overlap
        if (len(x)-fsize) % fshift == 0:
            post_pad = np.zeros((0,))
        else:
            post_pad = np.zeros((fshift - (len(x)-fsize) % fshift,))
        x = np.hstack((x, post_pad))
        M = (len(x)-fsize)//fshift + 1

            
    frame_starts = fshift * np.arange(M)    
    T=len(x)
    x = np.hstack((x, np.zeros(((M-1)*fshift + fsize - T,))))
    spec=np.zeros([M,fftsize//2+1]).astype('complex128')
    
    for m in xrange(M):
        frame = x[frame_starts[m]:frame_starts[m] + fsize]*awin
        temp  = scipy.fft(np.squeeze(frame),n=fftsize)
        spec[m]=temp[:fftsize//2+1]

    return spec


def istft(spec,fshift,swin,awin=None,fftsize=None,perfectrec=False):
    # iSTFT with a fixed frame shift

    if len(np.shape(spec)) != 2: # no multi-channel input
        raise ValueError('We only deal with single channel signals here')

    M, N = np.shape(spec)
    if N % 2 != 1:
        raise ValueError('We expect the spectrogram to only have non-negative frequencies')

    fsize = 2*(N-1)
    
    if awin is not None:
        # if an analysis window is passed, we assume the user wants us to normalize
        # the synthesis window accordingly. Otherwise, we trust the user.
        swin = synthwin(awin, fshift, swin=swin)

    if fftsize is None:
        fftsize = fsize

    if fftsize > len(swin):
        swin = np.hstack([swin,np.zeros((fftsize-len(swin),))])

    T= fshift * (M-1) + fsize
    signal=np.zeros(T)

    x_ran=np.arange(fsize)

    for s in xrange(M):
        iframe=np.real(scipy.ifft(np.concatenate((spec[s], spec[s][-2:0:-1].conjugate())),
                                  n=fftsize))

        iframe= iframe[0:fsize]
        signal[fshift*s+x_ran]+= iframe * np.squeeze(swin)

    if perfectrec is True:
        residual_size = fsize % fshift
        if residual_size == 0:
            pre_pad_length = fsize - fshift
        else:
            pre_pad_length = fsize - residual_size
        
        signal = signal[pre_pad_length:(fshift-fsize)]

    return signal


def get_consistency(S,fsize,fshift,awin,swin,perfectrec=False):
    # Compute the consistency
    # The user needs to make sure awin-swin are a perfect reconstruction pair, it is not enforced
    tmp = stft(istft(S,fshift,swin,perfectrec=perfectrec),fsize,fshift,awin,perfectrec=perfectrec)
    return 20 * np.log10(np.linalg.norm(S)/np.linalg.norm(tmp-S))

def extspec(S, L, Q):
    # Build an extended spectrograms to avoid requiring modulo in the computations
    T,Nreal = S.shape
    Np=Nreal+2*L
    Tp=T+2*(Q-1)
    ExtS = np.zeros((Tp,Np),dtype=S.dtype)
    ExtS[(Q-1):(Q-1+T),L:(Nreal+L)] = S
    ExtS[:,0:L] = np.conjugate(ExtS[:,(2*L):L:-1])
    ExtS[:,(Nreal+L):(Nreal+2*L)] = np.conjugate(ExtS[:,(Nreal+L-2):(Nreal-2):-1])
    ExtS[:(Q-1)]=np.atleast_2d(ExtS[Q-1])
    ExtS[(Q-1+T):] = np.atleast_2d(ExtS[Q-2+T])
    return ExtS
    

def create_weights(awin,swin,fshift,L):
    #  Compute the (L+1)xQxQ matrix of complex weights used by the LWS code
    T = len(awin)
    Q = np.int(T/fshift)
    interval = np.arange(L+1)
    expinterv= np.exp(-1j*2*np.pi*np.atleast_2d(interval).T*np.arange(T)/T)
    windowprod = np.zeros((T,Q))
    for q in range(Q):
        index=np.arange(T-q*fshift)
        windowprod[index,q] = awin[index] * swin[index+q*fshift]/T
    W = (expinterv.dot(windowprod)) * np.exp(-1j*2*np.pi*np.atleast_2d(interval).T*np.arange(Q)/Q)
    W[0,0] = W[0,0] - 1
    tmp = np.exp(1j*2*np.pi*np.atleast_2d(np.arange(Q)).T*np.arange(Q)/Q)
    W = W[:,np.newaxis] * tmp[np.newaxis,:]
    W = W.transpose((1,2,0)) # get back to row-major ordering, W is now Q*Q*(L+1)
    return W


def build_asymmetric_windows(awin_swin,fshift):
    # Compute the mirrored envelope used in Zhu et al.'s RTISI-LA.
    # Note that the input awin_swin should be the *product* of the analysis and
    # synthesis windows.
    T = len(awin_swin)
    Q = np.int(T/fshift)
    tmp = np.zeros((T,Q))
    tmp[:,0] = awin_swin
    for q in range(Q):
        index=np.arange(T-q*fshift)
        tmp[index,q] = awin_swin[q*fshift + index]
    
    win_ai = np.sum(tmp[:,1:],axis=1)[::-1]
    win_af = np.sum(tmp,axis=1)[::-1]
    if Q == 2:
        win_ai = awin_swin
    return win_ai, win_af


def get_thresholds(iterations,alpha,beta,gamma):
    # Compute the sparsity thresholds used to determine which T-F bins to update
    thresholds = alpha * np.exp(- beta * np.arange(iterations)**gamma)
    return thresholds


def batch_lws(S,W,thresholds):
    # Batch mode LWS phase reconstruction

    if S.dtype != np.complex128:
        S = S.astype(np.complex128)
    
    cdef int L = W.shape[2] - 1
    cdef int Q = W.shape[0]
    cdef int iterations = len(thresholds)
    if iterations == 0:
        return S
    cdef int T = S.shape[0]
    cdef int Nreal = S.shape[1]
    if Nreal % 2 == 0:
        raise ValueError('Please only include non-negative frequencies in the input spectrogram.')
    cdef int N = 2*(Nreal-1)
    
    cdef np.ndarray[np.double_t, ndim=3, mode="c"] Wr = np.ascontiguousarray(W.real)
    cdef np.ndarray[np.double_t, ndim=3, mode="c"] Wi = np.ascontiguousarray(W.imag)
    
    # Get a boolean mask specifying which weights to use or skip
    cdef double w_threshold = 1.0e-12
    cdef np.ndarray[int, ndim=3, mode="c"] Wflag = np.ascontiguousarray(np.abs(W) > w_threshold, dtype=np.dtype("i"))
    
    # Extend the spectrogram to avoid having to deal with values outside the boundaries
    ExtS  = extspec(S, L , Q)
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] ExtSr = np.ascontiguousarray(ExtS.real)
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] ExtSi = np.ascontiguousarray(ExtS.imag)
    # Store the amplitude spectrogram
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] AmpSpec = np.ascontiguousarray(np.abs(ExtS))
    cdef double mean_amp = np.mean(np.abs(S))

    # Perform the phase updates
    cdef double threshold
    for i in range(iterations):
        threshold = thresholds[i] * mean_amp
        if Q == 2:
            lwslib.LWSQ2(&ExtSr[0,0], &ExtSi[0,0], &Wr[0,0,0], &Wi[0,0,0], &Wflag[0,0,0], &AmpSpec[0,0], Nreal, T, L, threshold)
        elif Q == 4:
            lwslib.LWSQ4(&ExtSr[0,0], &ExtSi[0,0], &Wr[0,0,0], &Wi[0,0,0], &Wflag[0,0,0], &AmpSpec[0,0], Nreal, T, L, threshold)
        else:
            lwslib.LWSanyQ(&ExtSr[0,0], &ExtSi[0,0], &Wr[0,0,0], &Wi[0,0,0], &Wflag[0,0,0], &AmpSpec[0,0], Nreal, T, L, Q, threshold)
    
    # Extract the non-redundant part of the spectrogram
    S_out =  ExtSr[(Q-1):(Q-1+T),L:(Nreal+L)] + 1j * ExtSi[(Q-1):(Q-1+T),L:(Nreal+L)]
    
    return S_out


def nofuture_lws(S,W,thresholds):
    # Batch mode LWS phase reconstruction only considering past frames.
    # Typically only use for initialization.
    
    if S.dtype != np.complex128:
        S = S.astype(np.complex128)
    
    cdef int L = W.shape[2] - 1
    cdef int Q = W.shape[0]
    cdef int iterations = len(thresholds)
    if iterations == 0:
        return S
    cdef int T = S.shape[0]
    cdef int Nreal = S.shape[1]
    if Nreal % 2 == 0:
        raise ValueError('Please only include non-negative frequencies in the input spectrogram.')
    cdef int N = 2*(Nreal-1)
        
    cdef np.ndarray[np.double_t, ndim=3, mode="c"] Wr = np.ascontiguousarray(W.real)
    cdef np.ndarray[np.double_t, ndim=3, mode="c"] Wi = np.ascontiguousarray(W.imag)
    
    # Get a boolean mask specifying which weights to use or skip
    cdef double w_threshold = 1.0e-12
    cdef np.ndarray[int, ndim=3, mode="c"] Wflag = np.ascontiguousarray(np.abs(W) > w_threshold, dtype=np.dtype("i"))
    
    # Extend the spectrogram to avoid having to deal with values outside the boundaries
    ExtS  = extspec(S, L , Q)
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] ExtSr = np.ascontiguousarray(ExtS.real)
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] ExtSi = np.ascontiguousarray(ExtS.imag)
    # Store the amplitude spectrogram
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] AmpSpec = np.ascontiguousarray(np.abs(ExtS))
    cdef double mean_amp = np.mean(np.abs(S))
    
    # Perform the phase updates
    cdef double threshold
    for i in range(iterations):
        threshold = thresholds[i] * mean_amp
        if Q == 2:
            lwslib.NoFuture_LWSQ2(&ExtSr[0,0], &ExtSi[0,0], &Wr[0,0,0], &Wi[0,0,0], &Wflag[0,0,0], &AmpSpec[0,0], Nreal, T, L, threshold)
        elif Q == 4:
            lwslib.NoFuture_LWSQ4(&ExtSr[0,0], &ExtSi[0,0], &Wr[0,0,0], &Wi[0,0,0], &Wflag[0,0,0], &AmpSpec[0,0], Nreal, T, L, threshold)
        else:
            lwslib.NoFuture_LWSanyQ(&ExtSr[0,0], &ExtSi[0,0], &Wr[0,0,0], &Wi[0,0,0], &Wflag[0,0,0], &AmpSpec[0,0], Nreal, T, L,Q, threshold)
    
    # Extract the non-redundant part of the spectrogram
    S_out =  ExtSr[(Q-1):(Q-1+T),L:(Nreal+L)] + 1j * ExtSi[(Q-1):(Q-1+T),L:(Nreal+L)]
    
    return S_out


def online_lws(S,
               W, 
               W_ai, #W_asym_init, 
               W_af, #W_asym_full, 
               np.ndarray[np.double_t, ndim=1] thresholds,
               int LA):
    # Online mode LWS phase reconstruction
    
    if S.dtype != np.complex128:
        S = S.astype(np.complex128)

    cdef int L = W.shape[2] - 1
    cdef int Q = W.shape[0]
    cdef int iterations = len(thresholds)
    if iterations == 0:
        return S
    cdef int T = S.shape[0]
    cdef int Nreal = S.shape[1]
    if Nreal % 2 == 0:
        raise ValueError('Please only include non-negative frequencies in the input spectrogram.')
    cdef int N = 2*(Nreal-1)
    
    cdef np.ndarray[np.double_t, ndim=3, mode="c"] Wr = np.ascontiguousarray(W.real)
    cdef np.ndarray[np.double_t, ndim=3, mode="c"] Wi = np.ascontiguousarray(W.imag)
    cdef np.ndarray[np.double_t, ndim=3, mode="c"] Wr_ai = np.ascontiguousarray(W_ai.real)
    cdef np.ndarray[np.double_t, ndim=3, mode="c"] Wi_ai = np.ascontiguousarray(W_ai.imag)
    cdef np.ndarray[np.double_t, ndim=3, mode="c"] Wr_af = np.ascontiguousarray(W_af.real)
    cdef np.ndarray[np.double_t, ndim=3, mode="c"] Wi_af = np.ascontiguousarray(W_af.imag)
    
    # Get a boolean mask specifying which weights to use or skip
    cdef double w_threshold = 1.0e-12
    cdef np.ndarray[int, ndim=3, mode="c"] Wflag = np.ascontiguousarray(np.abs(W) > w_threshold, dtype=np.dtype("i"))
    cdef np.ndarray[int, ndim=3, mode="c"] Wflag_ai = np.ascontiguousarray(np.abs(W_ai) > w_threshold, dtype=np.dtype("i"))
    cdef np.ndarray[int, ndim=3, mode="c"] Wflag_af = np.ascontiguousarray(np.abs(W_af) > w_threshold, dtype=np.dtype("i"))
    
    # Extend the spectrogram to avoid having to deal with values outside the boundaries
    ExtS  = extspec(S, L , Q)
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] ExtSr = np.ascontiguousarray(ExtS.real)
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] ExtSi = np.ascontiguousarray(ExtS.imag)
    # Store the amplitude spectrogram
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] AmpSpec = np.ascontiguousarray(np.abs(ExtS))
    cdef double mean_amp = np.mean(np.abs(S))
    thresholds = np.ascontiguousarray(thresholds * mean_amp)

    cdef int update_type = 2
    # Perform the phase updates
    lwslib.TF_RTISI_LA(&ExtSr[0,0], &ExtSi[0,0], &Wr[0,0,0], &Wi[0,0,0],
                &Wr_ai[0,0,0], &Wi_ai[0,0,0], &Wr_af[0,0,0], &Wi_af[0,0,0],
                &Wflag[0,0,0], &Wflag_ai[0,0,0], &Wflag_af[0,0,0],
                &AmpSpec[0,0],
                iterations, LA, Nreal, T, L, Q, 
                &thresholds[0], update_type)

    # Extract the non-redundant part of the spectrogram
    S_out =  ExtSr[(Q-1):(Q-1+T),L:(Nreal+L)] + 1j * ExtSi[(Q-1):(Q-1+T),L:(Nreal+L)]
    
    return S_out


class lws(object):
    def __init__(self, awin_or_fsize, fshift, L = 5, swin = None, look_ahead = 3,
                 nofuture_iterations = 0, nofuture_alpha = 1, nofuture_beta = 0.1, nofuture_gamma = 1,
                 online_iterations = 0, online_alpha = 1, online_beta = 0.1, online_gamma = 1,
                 batch_iterations = 100, batch_alpha = 100, batch_beta = 0.1, batch_gamma = 1,
                 symmetric_win = True, mode= None, fftsize=None, perfectrec=True):
        if isinstance(awin_or_fsize, ( int, long ) ):
                # a frame size was passed in, build default perfect-reconstruction window
                awin = np.sqrt(hann(awin_or_fsize,symmetric=symmetric_win,use_offset=False))
                awin = np.sqrt( awin * synthwin(awin,fshift) )
        else:
            awin = awin_or_fsize
        if awin.ndim > 1:
            if (awin.ndim > 2) or (awin.shape[0]>1 and awin.shape>1):
                raise ValueError('The analysis window should be flat')
            else:
                awin = awin.flatten()
        
        if fftsize is None:
            fftsize = len(awin)

        if fftsize > len(awin):
            if (fftsize - len(awin)) % 2 != 0:
                raise ValueError('The zero-padding should add even length to the original window.')
            pad_length = (fftsize - len(awin)) // 2
            print('Zero-padding symmetrically around the original windows.\n'
                  'WARNING: for code simplicity, a consequence is that the first/last '
                  '{} samples of the signal will not be '.format(pad_length) + 
                  'in the perfect reconstruction region.')
            pad = np.zeros(pad_length)
            awin = np.hstack((pad,awin,pad))
            if swin is not None:
                swin = np.hstack((pad,swin,pad))
        
        self.awin = awin
        if swin is not None:
            print('Provided synthesis window is renormalized for perfect reconstruction.')
        self.swin = synthwin(awin,fshift,swin=swin)
        self.fshift = fshift
        self.fsize = len(awin)
        self.perfectrec = perfectrec
        self.L = L
        self.Q = np.int(self.fsize/self.fshift)
        self.W = create_weights(self.awin,self.swin,self.fshift,self.L)
        self.win_ai, self.win_af = build_asymmetric_windows(self.awin * self.swin, self.fshift)
        self.W_ai = create_weights(self.win_ai,self.swin,self.fshift,self.L)
        self.W_af = create_weights(self.win_af,self.swin,self.fshift,self.L)
        self.look_ahead = look_ahead
        
        if mode == 'speech':
            nofuture_iterations = 0
            online_iterations= 0
        elif mode == 'music':
            nofuture_iterations = 1
            online_iterations= 10

        self.batch_iterations = batch_iterations
        self.batch_alpha = batch_alpha
        self.batch_beta  = batch_beta
        self.batch_gamma = batch_gamma
        self.online_iterations = online_iterations
        self.online_alpha = online_alpha
        self.online_beta  = online_beta
        self.online_gamma = online_gamma
        self.nofuture_iterations = nofuture_iterations
        self.nofuture_alpha = nofuture_alpha
        self.nofuture_beta  = nofuture_beta
        self.nofuture_gamma = nofuture_gamma

        
        if (not np.allclose(awin, awin[::-1])):
            print('WARNING: It appears you are using an analysis window that is not symmetric.\n'
                  'The current code uses simplifications that rely on such symmetry, so the code may not behave properly.')


    def get_consistency(self,S):
        return get_consistency(S,self.fsize,self.fshift,self.awin,self.swin,perfectrec=self.perfectrec)


    def stft(self,S):
        return stft(S,self.fsize,self.fshift,self.awin,perfectrec=self.perfectrec)


    def istft(self,S):
        # we do not pass awin because swin was already renormalized for perfect reconstruction at initialization
        return istft(S,self.fshift,self.swin,perfectrec=self.perfectrec)


    def nofuture_lws(self,S,iterations=None,thresholds=None):
        if iterations is None:
            iterations = self.nofuture_iterations
        if thresholds is None:
            thresholds = get_thresholds(iterations,self.nofuture_alpha,self.nofuture_beta,self.nofuture_gamma)
        return nofuture_lws(S,self.W_ai,thresholds)


    def online_lws(self,S,iterations=None,thresholds=None):
        if iterations is None:
            iterations = self.online_iterations
        if thresholds is None:
            thresholds = get_thresholds(iterations,self.online_alpha,self.online_beta,self.online_gamma)
        return online_lws(S,self.W,self.W_ai,self.W_af,thresholds,self.look_ahead)


    def batch_lws(self,S,iterations=None,thresholds=None):
        if iterations is None:
            iterations = self.batch_iterations
        if thresholds is None:
            thresholds = get_thresholds(iterations,self.batch_alpha,self.batch_beta,self.batch_gamma)
        return batch_lws(S,self.W,thresholds)
        

    def run_lws(self,S):
        S0 = self.nofuture_lws(S)
        S1 = self.online_lws(S0)
        S2 = self.batch_lws(S1)
        return S2

