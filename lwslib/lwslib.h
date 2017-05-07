#ifndef LWSLIB_H_INCLUDED
#define LWSLIB_H_INCLUDED

#include <math.h>

void ExtendSpec(double * ExtSr, double * ExtSi, double * InSr, double * InSi, int Nreal, int M, int L, int Q);
void CopySpec(double * ExtSr, double * ExtSi, double * InSr, double * InSi, int Nreal, int M, int L, int Q);
void ComputeAmpSpec(double *Sr, double *Si, double * AmpSpec, int size);

void LWSQ2(double *Sr, double *Si, double *wr, double *wi, int *w_flag, double *AmpSpec, int Nreal, int M, int L, double threshold);
void LWSQ4(double *Sr, double *Si, double *wr, double *wi, int *w_flag, double *AmpSpec, int Nreal, int M, int L, double threshold);
void LWSanyQ(double *Sr, double *Si, double *wr, double *wi, int *w_flag, double *AmpSpec, int Nreal, int M, int L, int Q, double threshold);

void NoFuture_LWSQ2(double *Sr, double *Si, double *wr, double *wi, int *w_flag, double *AmpSpec, int Nreal, int M, int L, double threshold);
void NoFuture_LWSQ4(double *Sr, double *Si, double *wr, double *wi, int *w_flag, double *AmpSpec, int Nreal, int M, int L, double threshold);
void NoFuture_LWSanyQ(double *Sr, double *Si, double *wr, double *wi, int *w_flag, double *AmpSpec, int Nreal, int M, int L, int Q, double threshold);

void Asym_UpdatePhaseQ2(double *Sr, double *Si, double *wr, double *wi, int *w_flag, double *AmpSpec, int Nreal, int M, int M0, int L, double threshold, int update);
void Asym_UpdatePhaseQ4(double *Sr, double *Si, double *wr, double *wi, int *w_flag, double *AmpSpec, int Nreal, int M, int M0, int L, double threshold, int update);
void Asym_UpdatePhaseanyQ(double *Sr, double *Si, double *wr, double *wi, int *w_flag, double *AmpSpec, int Nreal, int M, int M0, int L, int Q, double threshold, int update);
void TF_RTISI_LA(double *Sr, double *Si, double *wr, double *wi, 
        double *wr_asym_init, double *wi_asym_init, double *wr_asym_full, double *wi_asym_full, int *w_flag, int *w_flag_ai, int *w_flag_af,  
        double *AmpSpec, int iter, int LA, int Nreal, int M, int L, int Q, double *ThresholdArray, int update);



#endif /* LWSLIB_H_INCLUDED */
