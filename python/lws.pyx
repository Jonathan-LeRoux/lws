cimport lwslib
cimport numpy as np
import numpy as np


def synthwin(awin,fshift,swin=None):
    # returns the normalized synthesis window for perfect reconstruction
    framesize = len(awin)
    Q = np.int(np.ceil(np.float(framesize) / np.float(fshift)))
    if swin is None:
        swin = awin
    twin = awin * swin                              
    tmp_framesize = Q*fshift
    
    w=np.hstack([twin,np.zeros((tmp_framesize-framesize,))])
    w=np.sum(np.reshape(w,(Q,fshift)),axis=0)
    w=np.tile(w,(1,Q))[0,:framesize]
    
    if min(w) <= 0:
        raise ValueError('The normalizer is not strictly positive')

    swin = swin / w

    return swin

def extspec(S, L, Q):
    Nreal,T = S.shape
    Np=Nreal+2*L
    Tp=T+2*(Q-1)
    ExtS = np.zeros((Np,Tp),dtype=S.dtype)
    ExtS[L:(Nreal+L),(Q-1):(Q-1+T)] = S
    ExtS[0:L] = np.conjugate(ExtS[(2*L):L:-1])
    ExtS[(Nreal+L):(Nreal+2*L)] = np.conjugate(ExtS[(Nreal+L-2):(Nreal-2):-1])
    ExtS[:,:(Q-1)]=np.atleast_2d(ExtS[:,Q-1]).T
    ExtS[:,(Q-1+T):] = np.atleast_2d(ExtS[:,Q-2+T]).T
    return ExtS
    
def create_weights(awin,swin,fshift,L):
    #  Compute the (L+1)xQxQ matrix of complex weights used by the LWS code
    Q = len(awin)/fshift
    # PLACEHOLDER:
    W = np.zeros((L+1,Q,Q),dtype = np.complex128)
    return W


def build_asymmetric_windows(awin_swin,fshift):
    # This computes the mirrored envelope used in Zhu et al.'s RTISI-LA.
    # Note that the input awin_swin should be the *product* of the analysis and
    # synthesis windows.
    Q = len(awin_swin)/fshift
    win_ai = awin_swin
    win_af = awin_swin
    return win_ai, win_af


def get_thresholds(iterations,alpha,beta,gamma):
    thresholds = alpha * np.exp(- beta * np.arange(iterations)**gamma)
    return thresholds


class lws(object):
    def __init__(self, awin, fshift, L = 2, swin = None, look_ahead = 3,
                 nofuture_iterations = 1, nofuture_alpha = 1, nofuture_beta = 0.1, nofuture_gamma = 1,
                 online_iterations = 10, online_alpha = 1, online_beta = 0.1, online_gamma = 1,
                 batch_iterations = 100, batch_alpha = 100, batch_beta = 0.1, batch_gamma = 1
                  ):
        if len(awin) % fshift != 0:
            raise ValueError('LWS requires that the window shift divides the window length.')

        self.awin = awin
        self.swin = synthwin(awin,fshift,swin=swin)
        self.fshift = fshift
        self.fsize = len(awin)
        self.L = L
        self.Q = self.fsize/self.fshift
        self.W = create_weights(self.awin,self.swin,self.fshift,self.L)
        self.W_ai, self.W_af = build_asymmetric_windows(self.awin * self.swin, self.fshift)
        self.batch_iterations = batch_iterations
        self.batch_alpha = batch_alpha
        self.batch_beta  = batch_beta
        self.batch_gamma = batch_gamma
        self.online_iterations = online_iterations
        self.online_alpha = online_alpha
        self.online_beta  = online_beta
        self.online_gamma = online_gamma
        self.look_ahead = look_ahead
        self.nofuture_iterations = nofuture_iterations
        self.nofuture_alpha = nofuture_alpha
        self.nofuture_beta  = nofuture_beta
        self.nofuture_gamma = nofuture_gamma


    def run_batch_lws(self,S,iterations=None,thresholds=None):
        if iterations is None:
            iterations = self.batch_iterations
        if thresholds is None:
            thresholds = get_thresholds(iterations,self.batch_alpha,self.batch_beta,self.batch_gamma)
        return batch_lws(S,self.W,thresholds)
        

    def run_online_lws(self,S,iterations=None,thresholds=None):
        if iterations is None:
            iterations = self.online_iterations
        if thresholds is None:
            thresholds = get_thresholds(iterations,self.online_alpha,self.online_beta,self.online_gamma)
        return online_lws(S,self.W,self.W_ai,self.W_af,thresholds,self.look_ahead)


    def run_nofuture_lws(self,S,iterations=None,thresholds=None):
        if iterations is None:
            iterations = self.nofuture_iterations
        if thresholds is None:
            thresholds = get_thresholds(iterations,self.nofuture_alpha,self.nofuture_beta,self.nofuture_gamma)
        return nofuture_lws(S,self.W_ai,thresholds)


    def run_lws(self,S):
        S0 = self.run_nofuture_lws(S)
        S1 = self.run_online_lws(S0)
        S2 = self.run_batch_lws(S1)
        return S2


def batch_lws(np.ndarray[np.complex128_t, ndim=2] S, 
              np.ndarray[np.complex128_t, ndim=3] W, 
              np.ndarray[np.double_t, ndim=1] thresholds):
        
    S = np.ascontiguousarray(S)
    W = np.ascontiguousarray(W)
    
    cdef int L = W.shape[0] - 1
    cdef int Q = W.shape[1]
    cdef int iterations = len(thresholds)
    cdef int Nreal = S.shape[0]
    cdef int T = S.shape[1]
    if Nreal % 2 == 0:
        raise ValueError('Please only include non-negative frequencies in the input spectrogram.')
    cdef int N = 2*(Nreal-1)
    
    cdef np.ndarray[np.double_t, ndim=3, mode="c"] Wr = np.ascontiguousarray(W.real)
    cdef np.ndarray[np.double_t, ndim=3, mode="c"] Wi = np.ascontiguousarray(W.imag)
    
    # Get a boolean mask specifying which weights to use or skip
    cdef double w_threshold = 1.0e-12
    cdef np.ndarray[int, ndim=3, mode="c"] Wflag = np.ascontiguousarray(np.abs(W) > w_threshold, dtype=np.dtype("i"))
    
    # Extend the spectrogram to avoid having to deal with values outside the boundaries
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] ExtS  = extspec(S, L , Q)
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
            lwslib.LWSanyQ(&ExtSr[0,0], &ExtSi[0,0], &Wr[0,0,0], &Wi[0,0,0], &Wflag[0,0,0], &AmpSpec[0,0], Nreal, T, L,Q, threshold)
    
    # Extract the non-redundant part of the spectrogram
    S_out =  ExtSr[L:(Nreal+L),(Q-1):(Q-1+T)] + 1j * ExtSi[L:(Nreal+L),(Q-1):(Q-1+T)]
    
    return S_out

def nofuture_lws(np.ndarray[np.complex128_t, ndim=2] S, 
        np.ndarray[np.complex128_t, ndim=3] W, 
        np.ndarray[np.double_t, ndim=1] thresholds):
        
    S = np.ascontiguousarray(S)
    W = np.ascontiguousarray(W)
    
    cdef int L = W.shape[0] - 1
    cdef int Q = W.shape[1]
    cdef int iterations = len(thresholds)
    cdef int Nreal = S.shape[0]
    cdef int T = S.shape[1]
    if Nreal % 2 == 0:
        raise ValueError('Please only include non-negative frequencies in the input spectrogram.')
    cdef int N = 2*(Nreal-1)
        
    cdef np.ndarray[np.double_t, ndim=3, mode="c"] Wr = np.ascontiguousarray(W.real)
    cdef np.ndarray[np.double_t, ndim=3, mode="c"] Wi = np.ascontiguousarray(W.imag)
    
    # Get a boolean mask specifying which weights to use or skip
    cdef double w_threshold = 1.0e-12
    cdef np.ndarray[int, ndim=3, mode="c"] Wflag = np.ascontiguousarray(np.abs(W) > w_threshold, dtype=np.dtype("i"))
    
    # Extend the spectrogram to avoid having to deal with values outside the boundaries
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] ExtS  = extspec(S, L , Q)
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
    S_out =  ExtSr[L:(Nreal+L),(Q-1):(Q-1+T)] + 1j * ExtSi[L:(Nreal+L),(Q-1):(Q-1+T)]
    
    return S_out


def online_lws(np.ndarray[np.complex128_t, ndim=2] S, 
        np.ndarray[np.complex128_t, ndim=3] W, 
        np.ndarray[np.complex128_t, ndim=3] W_asym_init, 
        np.ndarray[np.complex128_t, ndim=3] W_asym_full, 
        np.ndarray[np.double_t, ndim=1] thresholds,
	int LA):
        
        
    S          = np.ascontiguousarray(S)
    W          = np.ascontiguousarray(W)
    W_ai       = np.ascontiguousarray(W_asym_init)
    W_af       = np.ascontiguousarray(W_asym_full)
    
    cdef int L = W.shape[0] - 1
    cdef int Q = W.shape[1]
    cdef int iterations = len(thresholds)
    cdef int Nreal = S.shape[0]
    cdef int T = S.shape[1]
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
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] ExtS  = extspec(S, L , Q)
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
                &thresholds[0], update_type);

    # Extract the non-redundant part of the spectrogram
    S_out =  ExtSr[L:(Nreal+L),(Q-1):(Q-1+T)] + 1j * ExtSi[L:(Nreal+L),(Q-1):(Q-1+T)]
    
    return S_out
