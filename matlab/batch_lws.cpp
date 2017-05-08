/* LWS.CPP Construct consistent phase using Local Weighted Sums (LWS)
 * Syntax:    s_out = lws(s_in, weights, thresholds)
 *
 * Copyright (C) 2008-2017 Jonathan Le Roux
 * Apache 2.0  (http://www.apache.org/licenses/LICENSE-2.0) 
 */

#include "mex.h" /* Always include this */
#include "lwslib.h"
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
    #define THRESHOLDS prhs[2]
    if (nrhs < 3){
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
    
    
    if (nlhs > 0) {
        mxArray  *s_out;
        int T ,N, Nreal;
        double *pSr, *pSi, *pWr, *pWi, *pOr, *pOi, *pTr;
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
        
        pTr = mxGetPr(THRESHOLDS);
        
        int *pWflag;
        // Get a boolean mask specifying which weights to use or skip
        double w_threshold = 1.0e-12;
        pWflag = (int *) malloc(sizeof(int)*mxGetNumberOfElements(WEIGHTS));
        for(int n=0; n<mxGetNumberOfElements(WEIGHTS); n++){
            if (sqrt(pow(pWr[n], 2.)+pow(pWi[n], 2.)) > w_threshold){
                pWflag[n] = 1;
            } else {
                pWflag[n] = 0;
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
                
        // Perform the phase updates
        for (int j=0; j<iterations; j++){
            threshold = pTr[j] * mean_amp;
            //mexPrintf("%d: %f\n",j,threshold);
            if (Q==2){
                LWSQ2(ExtSr, ExtSi, pWr, pWi, pWflag, AmpSpec, Nreal, T, L, threshold);
            }
            else if (Q==4){
                LWSQ4(ExtSr, ExtSi, pWr, pWi, pWflag, AmpSpec, Nreal, T, L, threshold);
            }else {
                LWSanyQ(ExtSr, ExtSi, pWr, pWi, pWflag, AmpSpec, Nreal, T, L, Q, threshold);   
            }
        }
        
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
