cdef extern from "lwslib/lwslib.h":
    void LWSQ2(double *Sr, double *Si, double *wr, double *wi, int *w_flag, double *AmpSpec, int Nreal, int M, int L, double threshold)
    void LWSQ4(double *Sr, double *Si, double *wr, double *wi, int *w_flag, double *AmpSpec, int Nreal, int M, int L, double threshold)
    void LWSanyQ(double *Sr, double *Si, double *wr, double *wi, int *w_flag, double *AmpSpec, int Nreal, int M, int L, int Q, double threshold)
    void NoFuture_LWSQ2(double *Sr, double *Si, double *wr, double *wi, int *w_flag, double *AmpSpec, int Nreal, int M, int L, double threshold);
    void NoFuture_LWSQ4(double *Sr, double *Si, double *wr, double *wi, int *w_flag, double *AmpSpec, int Nreal, int M, int L, double threshold);
    void NoFuture_LWSanyQ(double *Sr, double *Si, double *wr, double *wi, int *w_flag, double *AmpSpec, int Nreal, int M, int L, int Q, double threshold);

    void TF_RTISI_LA(double *Sr, double *Si, double *wr, double *wi, 
        double *wr_asym_init, double *wi_asym_init, double *wr_asym_full, double *wi_asym_full, int *w_flag, int *w_flag_ai, int *w_flag_af,  
        double *AmpSpec, int iter, int LA, int Nreal, int M, int L, int Q, double *ThresholdArray, int update);


