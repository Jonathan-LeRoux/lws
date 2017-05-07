cimport lwslib
cimport numpy as np
import numpy as np


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
    

def lws(np.ndarray[np.complex128_t, ndim=2] S, 
        np.ndarray[np.complex128_t, ndim=3] W, 
        np.ndarray[np.double_t, ndim=1] thresholds):
        
        
    S          = np.ascontiguousarray(S)
    W          = np.ascontiguousarray(W)
    #thresholds = np.ascontiguousarray(thresholds)
    
    L = W.shape[0] - 1
    Q = W.shape[1]
    iterations = len(thresholds)
    Nreal = S.shape[0]
    T = S.shape[1]
    if Nreal % 2 == 0:
        raise ValueError('Please only include non-negative frequencies in the input spectrogram.')
    N = 2*(Nreal-1)
    
    cdef np.ndarray[np.double_t, ndim=3, mode="c"] Wr = np.ascontiguousarray(W.real)
    cdef np.ndarray[np.double_t, ndim=3, mode="c"] Wi = np.ascontiguousarray(W.imag)
    
    # Get a boolean mask specifying which weights to use or skip
    w_threshold = 1.0e-12
    cdef np.ndarray[int, ndim=3, mode="c"] Wflag = np.ascontiguousarray(np.abs(W) > w_threshold, dtype=np.dtype("i"))
    
    # Extend the spectrogram to avoid having to deal with values outside the boundaries
    ExtS = extspec(S, L , Q)
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] ExtSr = np.ascontiguousarray(ExtS.real)
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] ExtSi = np.ascontiguousarray(ExtS.imag)
    # Store the amplitude spectrogram
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] AmpSpec = np.ascontiguousarray(np.abs(ExtS))
    mean_amp = np.mean(np.abs(S))
    
    # Perform the phase updates
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
        
        
    S          = np.ascontiguousarray(S)
    W          = np.ascontiguousarray(W)
    #thresholds = np.ascontiguousarray(thresholds)
    
    L = W.shape[0] - 1
    Q = W.shape[1]
    iterations = len(thresholds)
    Nreal = S.shape[0]
    T = S.shape[1]
    if Nreal % 2 == 0:
        raise ValueError('Please only include non-negative frequencies in the input spectrogram.')
    N = 2*(Nreal-1)
    
    cdef np.ndarray[np.double_t, ndim=3, mode="c"] Wr = np.ascontiguousarray(W.real)
    cdef np.ndarray[np.double_t, ndim=3, mode="c"] Wi = np.ascontiguousarray(W.imag)
    
    # Get a boolean mask specifying which weights to use or skip
    w_threshold = 1.0e-12
    cdef np.ndarray[int, ndim=3, mode="c"] Wflag = np.ascontiguousarray(np.abs(W) > w_threshold, dtype=np.dtype("i"))
    
    # Extend the spectrogram to avoid having to deal with values outside the boundaries
    ExtS = extspec(S, L , Q)
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] ExtSr = np.ascontiguousarray(ExtS.real)
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] ExtSi = np.ascontiguousarray(ExtS.imag)
    # Store the amplitude spectrogram
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] AmpSpec = np.ascontiguousarray(np.abs(ExtS))
    mean_amp = np.mean(np.abs(S))
    
    # Perform the phase updates
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
    
    L = W.shape[0] - 1
    Q = W.shape[1]
    iterations = len(thresholds)
    Nreal = S.shape[0]
    T = S.shape[1]
    if Nreal % 2 == 0:
        raise ValueError('Please only include non-negative frequencies in the input spectrogram.')
    N = 2*(Nreal-1)
    
    cdef np.ndarray[np.double_t, ndim=3, mode="c"] Wr = np.ascontiguousarray(W.real)
    cdef np.ndarray[np.double_t, ndim=3, mode="c"] Wi = np.ascontiguousarray(W.imag)
    cdef np.ndarray[np.double_t, ndim=3, mode="c"] Wr_ai = np.ascontiguousarray(W_ai.real)
    cdef np.ndarray[np.double_t, ndim=3, mode="c"] Wi_ai = np.ascontiguousarray(W_ai.imag)
    cdef np.ndarray[np.double_t, ndim=3, mode="c"] Wr_af = np.ascontiguousarray(W_af.real)
    cdef np.ndarray[np.double_t, ndim=3, mode="c"] Wi_af = np.ascontiguousarray(W_af.imag)
    
    # Get a boolean mask specifying which weights to use or skip
    w_threshold = 1.0e-12
    cdef np.ndarray[int, ndim=3, mode="c"] Wflag = np.ascontiguousarray(np.abs(W) > w_threshold, dtype=np.dtype("i"))
    cdef np.ndarray[int, ndim=3, mode="c"] Wflag_ai = np.ascontiguousarray(np.abs(W_ai) > w_threshold, dtype=np.dtype("i"))
    cdef np.ndarray[int, ndim=3, mode="c"] Wflag_af = np.ascontiguousarray(np.abs(W_af) > w_threshold, dtype=np.dtype("i"))
    
    # Extend the spectrogram to avoid having to deal with values outside the boundaries
    ExtS = extspec(S, L , Q)
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] ExtSr = np.ascontiguousarray(ExtS.real)
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] ExtSi = np.ascontiguousarray(ExtS.imag)
    # Store the amplitude spectrogram
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] AmpSpec = np.ascontiguousarray(np.abs(ExtS))
    mean_amp = np.mean(np.abs(S))
    thresholds = np.ascontiguousarray(thresholds * mean_amp)

    update_type = 2
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
