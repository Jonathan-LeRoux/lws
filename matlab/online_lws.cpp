/* ONLINE_LWS.CPP Construct consistent phase using online LWS updates
 * Syntax:    s_out = online_lws(s_in, weights, thresholds)
 *
 * Copyright (C) 2008-2017 Jonathan Le Roux
 * Apache 2.0  (http://www.apache.org/licenses/LICENSE-2.0) 
 */

#include "mex.h" /* Always include this */
#include "lws_functions.h"
#include "matrix.h"

#define IS_REAL_1D_FULL_DOUBLE(P) (!mxIsComplex(P) && \
mxGetNumberOfDimensions(P) == 2 && !mxIsSparse(P) && mxIsDouble(P) && \
       (mxGetN(P)==1 || mxGetM(P)==1) )
#define IS_2D_FULL_DOUBLE(P) (mxGetNumberOfDimensions(P) == 2 && \
!mxIsSparse(P) && mxIsDouble(P))
#define IS_DOUBLE_SCALAR(P) (mxIsDouble(P) && !mxIsSparse(P) && mxGetNumberOfElements(P) == 1)



void mexFunction(int nlhs, mxArray *plhs[], /* Output variables */
        int nrhs, const mxArray *prhs[]) /* Input variables */ {
    
    /* Macros for the ouput and input arguments */
    #define S_IN prhs[0]
    #define WEIGHTS prhs[1]
    #define WEIGHTS_ASYM_INIT prhs[2]
    #define WEIGHTS_ASYM_FULL prhs[3]
    #define THRESHOLDS prhs[4]
    #define LA_IN prhs[5]
    if (nrhs < 6){
        mexPrintf("lws: not enought inputs\n");
        return;
    }      
    
    
    // Check spectrogram is correct type of input
    if (!IS_2D_FULL_DOUBLE(S_IN)) {
        mexPrintf("lws: spectrogram must be full 2-D double matrix.\n");
        return;
    }
    #ifdef DEBUG
    mexPrintf("lws: S_in is size %d x %d\n", mxGetM(S_IN), mxGetN(S_IN));
    #endif
            
    // Check weight matrix is correct type of input
    size_t K = mxGetNumberOfDimensions(WEIGHTS);
    const mwSize *Wdims = mxGetDimensions(WEIGHTS);
    if (K!= 3){
        mexPrintf("lws: weights should be 3-dimensional.\n");
    }
    int L = Wdims[0]-1;
    int Q = Wdims[1];

    #ifdef DEBUG
    mexPrintf("lws: weights is size %d x %d x %d\n", Wdims[0],Wdims[1],Wdims[2]);
    #endif
    // Check threshold list is correct type of input; also specifies the number of iterations
    if (!IS_REAL_1D_FULL_DOUBLE(THRESHOLDS)) {
        mexPrintf("lws: please provide a 1-D list of phase update thresholds.\n");
        return;
    }
    int iterations = (int) mxGetNumberOfElements(THRESHOLDS);
    #ifdef DEBUG
    mexPrintf("lws: will perform %d iterations with the specified thresholds.\n", iterations);
    #endif
    
    if (! IS_DOUBLE_SCALAR(LA_IN)){
        mexPrintf("Number of look-ahead frames is not a real scalar.\n");
        return;
    }
    int LA = (int) mxGetScalar(LA_IN);
    
    
    if (nlhs > 0) {
        mxArray  *s_out;
        int T ,N, Nreal;
        double *pSr, *pSi, *pWr, *pWi, *pOr, *pOi, *pTr;
        double *pWr_ai, *pWi_ai, *pWr_af, *pWi_af;
        double threshold;
        
        Nreal = (int) mxGetM(S_IN);
        T = (int) mxGetN(S_IN);
        
        if (Nreal%2 == 0){
            mexPrintf("Please only include non-negative frequencies in the input spectrogram.\n");
            return;
        }
        N = 2*(Nreal-1);
        
        // Set pointers to input data
        pSr = mxGetPr(S_IN);
        if (mxIsComplex(S_IN)) {
            pSi = mxGetPi(S_IN);
        }else {
            pSi = (double *) malloc(sizeof(double)*Nreal*T);
            for (int i=0; i<Nreal*T; i++) {
                pSi[i] = 0.;
            }
        }
        
        pWr = mxGetPr(WEIGHTS);
        pWi = mxGetPi(WEIGHTS);
        pWr_ai = mxGetPr(WEIGHTS_ASYM_INIT);
        pWi_ai = mxGetPi(WEIGHTS_ASYM_INIT);
        pWr_af = mxGetPr(WEIGHTS_ASYM_FULL);
        pWi_af = mxGetPi(WEIGHTS_ASYM_FULL);
        
        pTr = mxGetPr(THRESHOLDS);
        
        int *pWflag, *pWflag_ai, *pWflag_af;
        // Get a boolean mask specifying which weights to use or skip
        double w_threshold = 1.0e-12;
        pWflag    = (int *) malloc(sizeof(int)*mxGetNumberOfElements(WEIGHTS));
        pWflag_ai = (int *) malloc(sizeof(int)*mxGetNumberOfElements(WEIGHTS));
        pWflag_af = (int *) malloc(sizeof(int)*mxGetNumberOfElements(WEIGHTS));
        for(int n=0; n<mxGetNumberOfElements(WEIGHTS); n++){
            if (sqrt(pow(pWr[n], 2.)+pow(pWi[n], 2.)) > w_threshold){
                pWflag[n] = 1;
            } else {
                pWflag[n] = 0;
            }
            if (sqrt(pow(pWr_ai[n], 2.)+pow(pWi_ai[n], 2.)) > w_threshold){
                pWflag_ai[n] = 1;
            } else {
                pWflag_ai[n] = 0;
            }
            if (sqrt(pow(pWr_af[n], 2.)+pow(pWi_af[n], 2.)) > w_threshold){
                pWflag_af[n] = 1;
            } else {
                pWflag_af[n] = 0;
            }
        }

        // Extend the spectrogram to avoid having to deal with values outside the boundaries
        int Np=Nreal+2*L;
        double *ExtSr, *ExtSi;
        ExtSr = (double *) malloc(sizeof(double)*(T+2*(Q-1))*Np);
        ExtSi = (double *) malloc(sizeof(double)*(T+2*(Q-1))*Np);
        ExtendSpec(ExtSr, ExtSi, pSr, pSi, Nreal, T, L, Q);
        
        // Store the amplitude spectrogram
        double *AmpSpec;
        AmpSpec = (double *) malloc(sizeof(double)*(T+2*(Q-1))*Np);
        ComputeAmpSpec(ExtSr,ExtSi,AmpSpec,(T+2*(Q-1))*Np);
        double mean_amp = 0;
        for(int m=Q-1;m<(T+Q-1);m++){
            for(int n=0;n<Nreal;n++){
                mean_amp += AmpSpec[m*Np+n+L];
            }
        }
        mean_amp /= (T * Nreal);
                
        double *thresholds;
        thresholds = (double *) malloc(sizeof(double)*iterations);
        for (int j=0; j<iterations; j++){
            thresholds[j] = pTr[j] * mean_amp;
        }
        
        int update_type = 2;
        // Perform the phase updates
        TF_RTISI_LA(ExtSr, ExtSi, pWr, pWi, pWr_ai, pWi_ai, pWr_af, pWi_af, 
                pWflag, pWflag_ai, pWflag_af,
                AmpSpec, iterations, LA, Nreal, T, L, Q, thresholds, update_type);
        
        // Copy back the non-redundant part of the spectrogram into the output
        s_out = mxCreateDoubleMatrix(0, 0, mxCOMPLEX);
        mxSetM(s_out,Nreal);
        mxSetN(s_out,T);
        mxSetData(s_out, mxMalloc(sizeof(double)*T*Nreal));
        mxSetImagData(s_out, mxMalloc(sizeof(double)*T*Nreal));
        pOr = mxGetPr(s_out);
        pOi = mxGetPi(s_out);
        plhs[0] = s_out;
        CopySpec(ExtSr, ExtSi, pOr, pOi, Nreal, T, L, Q);
    }
    
    return;
}
