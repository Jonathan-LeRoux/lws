/* LWSLIB.CPP Core functions for constructing consistent phase 
 *                 using Local Weighted Sums (LWS)
 * 
 * Copyright (C) 2008-2018 Jonathan Le Roux
 * Apache 2.0  (http://www.apache.org/licenses/LICENSE-2.0) 
 */

#include "lwslib.h" 

/***************************************/
/*  Functions used by multiple files   */
/***************************************/

/* Build an extended spectrograms to avoid requiring modulo in the computations*/
void ExtendSpec(double * ExtSr, double * ExtSi, double * InSr, double * InSi, int Nreal, int M, int L, int Q) {
    int n, m, p;
    int Np=Nreal+2*L;
    
    for(m=0;m<(M+2*(Q-1));m++){
        p = m - Q + 1;
        if (p < 0){
            p = 0;
        } else if (p > M-1){
            p = M-1;
        }
        // copy negative frequencies from positive ones
        for(n=0;n<L;n++){
            ExtSr[m*Np+n] =  InSr[p*Nreal+L-n];
            ExtSi[m*Np+n] = -InSi[p*Nreal+L-n];
        }
        // copy the non-negative frequencies
        for(n=L;n<Nreal+L;n++){
            ExtSr[m*Np+n] =  InSr[p*Nreal+n-L];
            ExtSi[m*Np+n] =  InSi[p*Nreal+n-L];
        }
        // copy frequencies above Nyquist from those below
        for(n=Nreal+L;n<Nreal+2*L;n++){
            ExtSr[m*Np+n] =  ExtSr[m*Np+2*(Nreal+L-1)-n];
            ExtSi[m*Np+n] = -ExtSi[m*Np+2*(Nreal+L-1)-n];
        }
    }
    
    
}

/* Extract the original (non-redundant) part of the extended spectrogram */
void CopySpec(double * ExtSr, double * ExtSi, double * InSr, double * InSi, int Nreal, int M, int L, int Q) {
    int Np=Nreal+2*L;
    
    for(int m=0;m<M;m++){
        for(int n=0;n<Nreal;n++){
            InSr[m*Nreal+n]=ExtSr[(m+Q-1)*Np+n+L];
            InSi[m*Nreal+n]=ExtSi[(m+Q-1)*Np+n+L];
        }
    }
    
}

void ComputeAmpSpec(double *Sr, double *Si, double * AmpSpec, int size) {
    
    for(int n=0;n<size;n++){
        AmpSpec[n] = sqrt(pow(Sr[n], 2.)+pow(Si[n], 2.));
    }
    
}

/***************************************/
/* Functions used in lws.cpp           */
/***************************************/

// Using symmetries to reduce the number of calculations when Q=2
void LWSQ2(double *Sr, double *Si, double *wr, double *wi, int *w_flag, double *AmpSpec, int Nreal, int M, int L, double threshold) {
    int Q = 2;
    int n, m, k, r, mod;
    int Np=Nreal+2*L;
    double tempr, tempi, abstemp, absspec;
    double ar, ai, br, bi, cr, ci;
    
    int Naux=Nreal+L-1;
    
    for(m=Q-1;m<(M+Q-1);m++){
        for(n=L;n<Naux+1;n++){
            
            absspec = AmpSpec[m*Np+n];
            if(absspec>threshold){
                tempr=0.;
                tempi=0.;
                mod = ((n-L)%Q)*Q*(L+1);
                // contribution of n+/-k for center frame m
                for(k=1;k<L+1;k++){
                    if (w_flag[mod+k]) {
                        ar = wr[mod+k];
                        ai = wi[mod+k];
                        br = Sr[m*Np+n-k];
                        bi = Si[m*Np+n-k];
                        cr = Sr[m*Np+n+k];
                        ci = Si[m*Np+n+k];
                        tempr+= ar*(br+cr)-ai*(bi-ci);
                        tempi+= ar*(bi+ci)+ai*(br-cr);
                    }
                }
                // contribution of frames m+/-r
                for(r=1;r<Q;r++){
                    int u = r*(L+1);
                    int im = (m-r)*Np + n;
                    int ip = (m+r)*Np + n;
                    // for frequency n
                    if (w_flag[mod+u]) {
                        ar = wr[mod+u];
                        ai = wi[mod+u];
                        br = Sr[im];
                        bi = Si[im];
                        cr = Sr[ip];
                        ci = Si[ip];
                        tempr+= ar*(br+cr)-ai*(bi-ci);
                        tempi+= ar*(bi+ci)+ai*(br-cr);
                    }
                    // for frequencies n+/-k
                    for(k=1;k<L+1;k++){
                        if (w_flag[mod+u+k]) {
                            ar = wr[mod+u+k];
                            ai = wi[mod+u+k];
                            br = Sr[im-k] + Sr[ip+k];
                            bi = Si[im-k] + Si[ip+k];
                            cr = Sr[ip-k] + Sr[im+k];
                            ci = Si[ip-k] + Si[im+k];
                            tempr+= ar*(br+cr)-ai*(bi-ci);
                            tempi+= ar*(bi+ci)+ai*(br-cr);
                        }
                    }
                }
                
                abstemp = sqrt(pow(tempr, 2.)+pow(tempi, 2.));
                if((abstemp>0)){
                    // update the phase
                    Sr[m*Np+n]= tempr*absspec/abstemp;
                    Si[m*Np+n]= tempi*absspec/abstemp;
                    // propagate changes in the extended regions
                    if((n>=L+1)&&(n<2*L+1)){
                        Sr[m*Np+2*L-n]= Sr[m*Np+n];
                        Si[m*Np+2*L-n]= -Si[m*Np+n];
                    }else if((n>=Nreal-1)&&(n<Naux)){
                        Sr[m*Np+2*Naux-n] = Sr[m*Np+n];
                        Si[m*Np+2*Naux-n] = -Si[m*Np+n];
                    }
                }
            }
        }
    }
}

// Using symmetries to reduce the number of calculations when Q=4
void LWSQ4(double *Sr, double *Si, double *wr, double *wi, int *w_flag, double *AmpSpec, int Nreal, int M, int L, double threshold) {
    int Q = 4;
    int n, m, k, r, mod;
    int Np=Nreal+2*L;
    double tempr, tempi, abstemp, absspec;
    double ar, ai, br, bi, cr, ci;
    
    int Naux=Nreal+L-1;
    
    for(m=Q-1;m<(M+Q-1);m++){
        for(n=L;n<Naux+1;n++){
            
            absspec = AmpSpec[m*Np+n];
            if(absspec>threshold){
                
                tempr=0.;
                tempi=0.;
                mod = ((n-L)%Q)*Q*(L+1);
                // contribution of n+/-k for center frame m
                for(k=1;k<L+1;k++){
                    if (w_flag[mod+k]) {
                        ar = wr[mod+k];
                        ai = wi[mod+k];
                        br = Sr[m*Np+n-k];
                        bi = Si[m*Np+n-k];
                        cr = Sr[m*Np+n+k];
                        ci = Si[m*Np+n+k];
                        tempr+= ar*(br+cr)-ai*(bi-ci);
                        tempi+= ar*(bi+ci)+ai*(br-cr);
                    }
                }
                // contribution of frames m+/-r
                // different simplifications can be made depending on parity
                if((n-L)%2==1){
                    for(r=1;r<Q;r+=2){
                        // for frequency n
                        if (w_flag[mod+r*(L+1)]) {
                            ar = wr[mod+r*(L+1)];
                            ai = wi[mod+r*(L+1)];
                            br = Sr[(m-r)*Np+n];
                            bi = Si[(m-r)*Np+n];
                            cr = Sr[(m+r)*Np+n];
                            ci = Si[(m+r)*Np+n];
                            tempr+= ar*(br+cr)-ai*(bi-ci);
                            tempi+= ar*(bi+ci)+ai*(br-cr);
                        }
                        // for frequencies n+/-k
                        for(k=1;k<L+1;k++){
                            if (w_flag[mod+r*(L+1)+k]) {
                                ar = wr[mod+r*(L+1)+k];
                                ai = wi[mod+r*(L+1)+k];
                                br = Sr[(m-r)*Np+n-k] - Sr[(m+r)*Np+n+k];
                                bi = Si[(m-r)*Np+n-k] - Si[(m+r)*Np+n+k];
                                cr = Sr[(m+r)*Np+n-k] - Sr[(m-r)*Np+n+k];
                                ci = Si[(m+r)*Np+n-k] - Si[(m-r)*Np+n+k];
                                tempr+= ar*(br+cr)-ai*(bi-ci);
                                tempi+= ar*(bi+ci)+ai*(br-cr);
                            }
                        }
                    }
                    r=2;
                    if (w_flag[mod+r*(L+1)]) {
                        ar = wr[mod+r*(L+1)];
                        ai = wi[mod+r*(L+1)];
                        br = Sr[(m-r)*Np+n];
                        bi = Si[(m-r)*Np+n];
                        cr = Sr[(m+r)*Np+n];
                        ci = Si[(m+r)*Np+n];
                        tempr+= ar*(br+cr)-ai*(bi-ci);
                        tempi+= ar*(bi+ci)+ai*(br-cr);
                    }
                    for(k=1;k<L+1;k++){
                        if (w_flag[mod+r*(L+1)+k]) {
                            ar = wr[mod+r*(L+1)+k];
                            ai = wi[mod+r*(L+1)+k];
                            br = Sr[(m-r)*Np+n-k] + Sr[(m+r)*Np+n+k];
                            bi = Si[(m-r)*Np+n-k] + Si[(m+r)*Np+n+k];
                            cr = Sr[(m+r)*Np+n-k] + Sr[(m-r)*Np+n+k];
                            ci = Si[(m+r)*Np+n-k] + Si[(m-r)*Np+n+k];
                            tempr+= ar*(br+cr)-ai*(bi-ci);
                            tempi+= ar*(bi+ci)+ai*(br-cr);
                        }
                    }
                }else{
                    for(r=1;r<Q;r++){
                        if (w_flag[mod+r*(L+1)]) {
                            ar = wr[mod+r*(L+1)];
                            ai = wi[mod+r*(L+1)];
                            br = Sr[(m-r)*Np+n];
                            bi = Si[(m-r)*Np+n];
                            cr = Sr[(m+r)*Np+n];
                            ci = Si[(m+r)*Np+n];
                            tempr+= ar*(br+cr)-ai*(bi-ci);
                            tempi+= ar*(bi+ci)+ai*(br-cr);
                        }
                        for(k=1;k<L+1;k++){
                            if (w_flag[mod+r*(L+1)+k]) {
                                ar = wr[mod+r*(L+1)+k];
                                ai = wi[mod+r*(L+1)+k];
                                br = Sr[(m-r)*Np+n-k] + Sr[(m+r)*Np+n+k];
                                bi = Si[(m-r)*Np+n-k] + Si[(m+r)*Np+n+k];
                                cr = Sr[(m+r)*Np+n-k] + Sr[(m-r)*Np+n+k];
                                ci = Si[(m+r)*Np+n-k] + Si[(m-r)*Np+n+k];
                                tempr+= ar*(br+cr)-ai*(bi-ci);
                                tempi+= ar*(bi+ci)+ai*(br-cr);
                            }
                        }
                    }
                }
                
                abstemp = sqrt(pow(tempr, 2.)+pow(tempi, 2.));
                if((abstemp>0)){
                    // update the phase
                    Sr[m*Np+n]= tempr*absspec/abstemp;
                    Si[m*Np+n]= tempi*absspec/abstemp;
                    // propagate changes in the extended regions
                    if((n>=L+1)&&(n<2*L+1)){
                        Sr[m*Np+2*L-n]= Sr[m*Np+n];
                        Si[m*Np+2*L-n]= -Si[m*Np+n];
                    }else if((n>=Nreal-1)&&(n<Naux)){
                        Sr[m*Np+2*Naux-n] = Sr[m*Np+n];
                        Si[m*Np+2*Naux-n] = -Si[m*Np+n];
                    }
                }
            }
        }
    }
}

// General case, slightly slower for Q=2 or Q=4
void LWSanyQ(double *Sr, double *Si, double *wr, double *wi, int *w_flag, double *AmpSpec, int Nreal, int M, int L, int Q, double threshold) {
    int n, m, k, r, mod, modneg;
    int Np=Nreal+2*L;
    double tempr, tempi, abstemp, absspec;
    double ar, ai, br, bi, cr, ci;
    
    int Naux=Nreal+L-1;
    
    for(m=Q-1;m<(M+Q-1);m++){
        for(n=L;n<Naux+1;n++){
            
            // only update bins whose magnitude is above threshold
            absspec = AmpSpec[m*Np+n];
            if(absspec>threshold){
                tempr=0.;
                tempi=0.;
                mod = ((n-L)%Q)*Q*(L+1);
                modneg = ((Q-((n-L)%Q))%Q)*Q*(L+1);
                // contribution of n+/-k for center frame m
                for(k=1;k<L+1;k++){
                    if (w_flag[mod+k]) {
                        ar = wr[mod+k];
                        ai = wi[mod+k];
                        br = Sr[m*Np+n-k];
                        bi = Si[m*Np+n-k];
                        cr = Sr[m*Np+n+k];
                        ci = Si[m*Np+n+k];
                        tempr+= ar*(br+cr)-ai*(bi-ci);
                        tempi+= ar*(bi+ci)+ai*(br-cr);
                    }
                }
                
                // contribution of frames m+/-r
                for(r=1;r<Q;r++){
                    int u = r*(L+1);
                    int im = (m-r)*Np + n;
                    int ip = (m+r)*Np + n;
                    // for frequency n
                    if (w_flag[mod+u]) {
                        ar = wr[mod+u];
                        ai = wi[mod+u];
                        br = Sr[im];
                        bi = Si[im];
                        cr = Sr[ip];
                        ci = Si[ip];
                        tempr+= ar*(br+cr)-ai*(bi-ci);
                        tempi+= ar*(bi+ci)+ai*(br-cr);
                    }
                    // for frequencies n+/-k
                    for(k=1;k<L+1;k++){
                        if (w_flag[mod+u+k]) {
                            ar = wr[mod+u+k];
                            ai = wi[mod+u+k];
                            br = Sr[im-k];
                            bi = Si[im-k];
                            cr = Sr[ip-k];
                            ci = Si[ip-k];
                            tempr+= ar*(br+cr)-ai*(bi-ci);
                            tempi+= ar*(bi+ci)+ai*(br-cr);
                        }
                        if (w_flag[modneg+u+k]) {
                            ar = wr[modneg+u+k];
                            ai = wi[modneg+u+k];
                            br = Sr[ip+k];
                            bi = Si[ip+k];
                            cr = Sr[im+k];
                            ci = Si[im+k];
                            tempr+= ar*(br+cr)-ai*(bi-ci);
                            tempi+= ar*(bi+ci)+ai*(br-cr);
                        }
                    }
                }
                
                abstemp = sqrt(pow(tempr, 2.)+pow(tempi, 2.));
                if((abstemp>0)){
                    // update the phase
                    Sr[m*Np+n]= tempr*absspec/abstemp;
                    Si[m*Np+n]= tempi*absspec/abstemp;
                    // propagate changes in the extended regions
                    if((n>=L+1)&&(n<2*L+1)){
                        Sr[m*Np+2*L-n]= Sr[m*Np+n];
                        Si[m*Np+2*L-n]= -Si[m*Np+n];
                    }else if((n>=Nreal-1)&&(n<Naux)){
                        Sr[m*Np+2*Naux-n] = Sr[m*Np+n];
                        Si[m*Np+2*Naux-n] = -Si[m*Np+n];
                    }
                }
            }
        }
    }
}


/***************************************/
/* Functions used by nofuture_lws.cpp  */
/***************************************/

void NoFuture_LWSQ2(double *Sr, double *Si, double *wr, double *wi, int *w_flag, double *AmpSpec, int Nreal, int M, int L, double threshold) {
    int Q = 2;
    int n, m, k, r, mod;
    int Np=Nreal+2*L;
    double tempr, tempi, abstemp, absspec;
    double ar, ai, br, bi, cr, ci;
    
    int Naux=Nreal+L-1;
    
    for(m=Q-1;m<(M+Q-1);m++){
        for(n=L;n<Naux+1;n++){
            
            absspec = AmpSpec[m*Np+n];
            if(absspec>threshold){
                tempr=0.;
                tempi=0.;
                mod = ((n-L)%Q)*Q*(L+1);
                // contribution of frames m-r
                for(r=1;r<Q;r++){
                    int u = r*(L+1);
                    int im = (m-r)*Np + n;
                    // for frequency n
                    if (w_flag[mod+u]) {
                        ar = wr[mod+u];
                        ai = wi[mod+u];
                        br = Sr[im];
                        bi = Si[im];
                        tempr+= ar*(br)-ai*(bi);
                        tempi+= ar*(bi)+ai*(br);
                    }
                    // for frequencies n+/-k
                    for(k=1;k<L+1;k++){
                        if (w_flag[mod+u+k]) {
                            ar = wr[mod+u+k];
                            ai = wi[mod+u+k];
                            br = Sr[im-k];
                            bi = Si[im-k];
                            cr = Sr[im+k];
                            ci = Si[im+k];
                            tempr+= ar*(br+cr)-ai*(bi-ci);
                            tempi+= ar*(bi+ci)+ai*(br-cr);
                        }
                    }
                }
                
                abstemp = sqrt(pow(tempr, 2.)+pow(tempi, 2.));
                if((abstemp>0)){
                    // update the phase
                    Sr[m*Np+n]= tempr*absspec/abstemp;
                    Si[m*Np+n]= tempi*absspec/abstemp;
                    // propagate changes in the extended regions
                    if((n>=L+1)&&(n<2*L+1)){
                        Sr[m*Np+2*L-n]= Sr[m*Np+n];
                        Si[m*Np+2*L-n]= -Si[m*Np+n];
                    }else if((n>=Nreal-1)&&(n<Naux)){
                        Sr[m*Np+2*Naux-n] = Sr[m*Np+n];
                        Si[m*Np+2*Naux-n] = -Si[m*Np+n];
                    }
                }
            }
        }
    }
}

// Using symmetries to reduce the number of calculations when Q=4
void NoFuture_LWSQ4(double *Sr, double *Si, double *wr, double *wi, int *w_flag, double *AmpSpec, int Nreal, int M, int L, double threshold) {
    int Q = 4;
    int n, m, k, r, mod;
    int Np=Nreal+2*L;
    double tempr, tempi, abstemp, absspec;
    double ar, ai, br, bi, cr, ci;
    
    int Naux=Nreal+L-1;
    
    for(m=Q-1;m<(M+Q-1);m++){
        for(n=L;n<Naux+1;n++){
            
            absspec = AmpSpec[m*Np+n];
            if(absspec>threshold){
                
                tempr=0.;
                tempi=0.;
                mod = ((n-L)%Q)*Q*(L+1);
                
                for(r=Q-1;r>0;r--){
                    int u = r*(L+1);
                    int im = (m-r)*Np + n;
                    // contribution of frames m-r for frequencies n+/-k
                    // different simplifications can be made depending on parity
                    if((n-L)%2==1 && r%2==1){
                        for(k=1;k<L+1;k++){
                            if (w_flag[mod+u+k]) {
                                ar = wr[mod+u+k];
                                ai = wi[mod+u+k];
                                br = Sr[im+n-k];
                                bi = Si[im+n-k];
                                cr = - Sr[im+n+k];
                                ci = - Si[im+n+k];
                                tempr+= ar*(br+cr)-ai*(bi-ci);
                                tempi+= ar*(bi+ci)+ai*(br-cr);
                            }
                        }
                    } else {
                        for(k=1;k<L+1;k++){
                            if (w_flag[mod+u+k]) {
                                ar = wr[mod+u+k];
                                ai = wi[mod+u+k];
                                br = Sr[im+n-k];
                                bi = Si[im+n-k];
                                cr = Sr[im+n+k];
                                ci = Si[im+n+k];
                                tempr+= ar*(br+cr)-ai*(bi-ci);
                                tempi+= ar*(bi+ci)+ai*(br-cr);
                            }
                        }
                    }
                    // contribution of frames m-Q+1 to m-1 for frequency n
                    if (w_flag[mod+u]) {
                        ar = wr[mod+u];
                        ai = wi[mod+u];
                        br = Sr[im+n];
                        bi = Si[im+n];
                        tempr+= ar*(br)-ai*(bi);
                        tempi+= ar*(bi)+ai*(br);
                    }
                }

                abstemp = sqrt(pow(tempr, 2.)+pow(tempi, 2.));
                if((abstemp>0)){
                    // update the phase
                    Sr[m*Np+n]= tempr*absspec/abstemp;
                    Si[m*Np+n]= tempi*absspec/abstemp;
                    // propagate changes in the extended regions
                    if((n>=L+1)&&(n<2*L+1)){
                        Sr[m*Np+2*L-n]= Sr[m*Np+n];
                        Si[m*Np+2*L-n]= -Si[m*Np+n];
                    }else if((n>=Nreal-1)&&(n<Naux)){
                        Sr[m*Np+2*Naux-n] = Sr[m*Np+n];
                        Si[m*Np+2*Naux-n] = -Si[m*Np+n];
                    }
                }
            }
        }
    }
}

// General case, slightly slower for Q=2 or Q=4
void NoFuture_LWSanyQ(double *Sr, double *Si, double *wr, double *wi, int *w_flag, double *AmpSpec, int Nreal, int M, int L, int Q, double threshold) {
    int n, m, k, r, mod, modneg;
    int Np=Nreal+2*L;
    double tempr, tempi, abstemp, absspec;
    double ar, ai, br, bi, cr, ci;
    
    int Naux=Nreal+L-1;
    
    for(m=Q-1;m<(M+Q-1);m++){
        for(n=L;n<Naux+1;n++){
            
            // only update bins whose magnitude is above threshold
            absspec = AmpSpec[m*Np+n];
            if(absspec>threshold){
                tempr=0.;
                tempi=0.;
                mod = ((n-L)%Q)*Q*(L+1);
                modneg = ((Q-((n-L)%Q))%Q)*Q*(L+1);
                
                // contribution of frames m-r
                for(r=1;r<Q;r++){
                    int u = r*(L+1);
                    int im = (m-r)*Np + n;
                    // for frequency n
                    if (w_flag[mod+u]) {
                        ar = wr[mod+u];
                        ai = wi[mod+u];
                        br = Sr[im];
                        bi = Si[im];
                        tempr+= ar*(br)-ai*(bi);
                        tempi+= ar*(bi)+ai*(br);
                    }
                    // for frequencies n+/-k
                    for(k=1;k<L+1;k++){
                        if (w_flag[mod+u+k]) {
                            ar = wr[mod+u+k];
                            ai = wi[mod+u+k];
                            br = Sr[im-k];
                            bi = Si[im-k];
                            tempr+= ar*(br)-ai*(bi);
                            tempi+= ar*(bi)+ai*(br);
            }
                if (w_flag[modneg+u+k]) {
                            ar = wr[modneg+u+k];
                            ai = wi[modneg+u+k];
                            cr = Sr[im+k];
                            ci = Si[im+k];
                            tempr+= ar*(cr)+ai*(ci);
                            tempi+= ar*(ci)-ai*(cr);
                        }
                    }
                }
                
                abstemp = sqrt(pow(tempr, 2.)+pow(tempi, 2.));
                if((abstemp>0)){
                    // update the phase
                    Sr[m*Np+n]= tempr*absspec/abstemp;
                    Si[m*Np+n]= tempi*absspec/abstemp;
                    // propagate changes in the extended regions
                    if((n>=L+1)&&(n<2*L+1)){
                        Sr[m*Np+2*L-n]= Sr[m*Np+n];
                        Si[m*Np+2*L-n]= -Si[m*Np+n];
                    }else if((n>=Nreal-1)&&(n<Naux)){
                        Sr[m*Np+2*Naux-n] = Sr[m*Np+n];
                        Si[m*Np+2*Naux-n] = -Si[m*Np+n];
                    }
                }
            }
        }
    }
}


/***************************************/
/* Functions used by online_lws.cpp    */
/***************************************/

// M gives the number of frames to update, and M0 the number of frames on the right 
// to use in the computations. 
// Frames 0 to M-1 are updated using data from frames -Q+1 to M0.
//

void Asym_UpdatePhaseQ2(double *Sr, double *Si, double *wr, double *wi, int *w_flag, double *AmpSpec, int Nreal, int M, int M0, int L, double threshold, int update)
{
    int Q = 2;
    int n,m,k,r,mod;
    int Np=Nreal+2*L;
    double tempr, tempi, abstemp, absspec;
    double ar,ai,br,bi,cr,ci;

    int Naux=Nreal+L-1;
    int rframe=0; // indicates the last frame to use on the right
    int cframe=0;
    
    for(m=Q-1;m<(M+Q-1);m++){
        
        rframe= M0+Q-m-1;
        if(rframe>Q){
            rframe=Q;
        }
        cframe=1;
        if(rframe<1){
            cframe=0;
            rframe=1;
        }
        
        for(n=L;n<Naux+1;n++){
            
            absspec = AmpSpec[m*Np+n];
            if(absspec>threshold){
                tempr=0.;
                tempi=0.;
                mod = ((n-L)%Q)*Q*(L+1);
                
                // contribution of center frame m
                if(cframe){
                    if(update==1){
                        tempr+=Sr[m*Np+n]/Q;
                        tempi+=Si[m*Np+n]/Q;
                    }
                    for(k=1;k<L+1;k++){
                        if (w_flag[mod+k]) {
                            ar = wr[mod+k];
                            ai = wi[mod+k];
                            br = Sr[m*Np+n-k];
                            bi = Si[m*Np+n-k];
                            cr = Sr[m*Np+n+k];
                            ci = Si[m*Np+n+k];
                            tempr+= ar*(br+cr)-ai*(bi-ci);
                            tempi+= ar*(bi+ci)+ai*(br-cr);
                        }
                    }
                }
                
                for(r=1;r<rframe;r++){// add frames  at m+/-r (i.e., on left and right)
                    int u = r*(L+1);
                    int im = (m-r)*Np + n;
                    int ip = (m+r)*Np + n;
                    // for frequency n
                    if (w_flag[mod+u]) {
                        ar = wr[mod+u];
                        ai = wi[mod+u];
                        br = Sr[im];
                        bi = Si[im];
                        cr = Sr[ip];
                        ci = Si[ip];
                        tempr+= ar*(br+cr)-ai*(bi-ci);
                        tempi+= ar*(bi+ci)+ai*(br-cr);  
                    }
                    for(k=1;k<L+1;k++){
                        if (w_flag[mod+u+k]) {
                            ar = wr[mod+u+k];
                            ai = wi[mod+u+k];
                            br = Sr[im-k] + Sr[ip+k];
                            bi = Si[im-k] + Si[ip+k];
                            cr = Sr[ip-k] + Sr[im+k];
                            ci = Si[ip-k] + Si[im+k];
                            tempr+= ar*(br+cr)-ai*(bi-ci);
                            tempi+= ar*(bi+ci)+ai*(br-cr);
                        }
                    }
                }
                
                for(r=rframe;r<Q;r++){ // for those only the left frames should be added
                    int u = r*(L+1);
                    int im = (m-r)*Np + n;
                    if (w_flag[mod+u]) {
                        ar = wr[mod+u];
                        ai = wi[mod+u];
                        br = Sr[im];
                        bi = Si[im];
                        tempr+= ar*br-ai*bi;
                        tempi+= ar*bi+ai*br;
                    }    
                    for(k=1;k<L+1;k++){
                        if (w_flag[mod+u+k]) {
                            ar = wr[mod+u+k];
                            ai = wi[mod+u+k];
                            br = Sr[im-k];
                            bi = Si[im-k];
                            cr = Sr[im+k];
                            ci = Si[im+k];
                            tempr+= ar*(br+cr)-ai*(bi-ci);
                            tempi+= ar*(bi+ci)+ai*(br-cr);
                        }
                    }
                }
                
                abstemp = sqrt(pow(tempr, 2.)+pow(tempi, 2.));
                if((abstemp>0)){
                    // update the phase
                    Sr[m*Np+n]= tempr*absspec/abstemp;
                    Si[m*Np+n]= tempi*absspec/abstemp;
                    // propagate changes in the extended regions
                    if((n>=L+1)&&(n<2*L+1)){
                        Sr[m*Np+2*L-n]= Sr[m*Np+n];
                        Si[m*Np+2*L-n]= -Si[m*Np+n];
                    }else if((n>=Nreal-1)&&(n<Naux)){
                        Sr[m*Np+2*Naux-n] = Sr[m*Np+n];
                        Si[m*Np+2*Naux-n] = -Si[m*Np+n];
                    }
                }
            }
        }
    }
}

// Same as above for Q=4
void Asym_UpdatePhaseQ4(double *Sr, double *Si, double *wr, double *wi, int *w_flag, double *AmpSpec, int Nreal, int M, int M0, int L, double threshold, int update)
{
    int Q = 4;
    int n,m,k,r,mod;
    int Np=Nreal+2*L;
    double tempr, tempi, abstemp, absspec;
    double ar,ai,br,bi,cr,ci;

    int Naux=Nreal+L-1;
    int rframe=0; // indicates the last frame to use on the right
    int cframe=0;

    
    for(m=Q-1;m<(M+Q-1);m++){
        
        rframe= M0+Q-m-1;
        if(rframe>Q){
            rframe=Q;
        }
        cframe=1;
        if(rframe<1){
            cframe=0;
            rframe=1;
        }
        
        for(n=L;n<Naux+1;n++){
            
            absspec = AmpSpec[m*Np+n];
            if(absspec>threshold){
                tempr=0.;
                tempi=0.;
                mod = ((n-L)%Q)*Q*(L+1);
        // contribution of center frame m
                if(cframe){
                    if(update==1){
                        tempr+=Sr[m*Np+n]/Q;
                        tempi+=Si[m*Np+n]/Q;
                    }
                    for(k=1;k<L+1;k++){
                        if (w_flag[mod+k]) {
                            ar = wr[mod+k];
                            ai = wi[mod+k];
                            br = Sr[m*Np+n-k];
                            bi = Si[m*Np+n-k];
                            cr = Sr[m*Np+n+k];
                            ci = Si[m*Np+n+k];
                            tempr+= ar*(br+cr)-ai*(bi-ci);
                            tempi+= ar*(bi+ci)+ai*(br-cr);
                        }
                    }
                }
                if((n-L)%2==1){
                    for(r=1;r<Q;r+=2){
                        int u = r*(L+1);
                        int im = (m-r)*Np + n;
                        int ip = (m+r)*Np + n;
                        if(r<rframe){
                            if (w_flag[mod+u]) {
                                ar = wr[mod+u];
                                ai = wi[mod+u];
                                br = Sr[im];
                                bi = Si[im];
                                cr = Sr[ip];
                                ci = Si[ip];
                                tempr+= ar*(br+cr)-ai*(bi-ci);
                                tempi+= ar*(bi+ci)+ai*(br-cr);
                            }
                            for(k=1;k<L+1;k++){
                                if (w_flag[mod+u+k]) {
                                    ar = wr[mod+u+k];
                                    ai = wi[mod+u+k];
                                    br = Sr[im-k] - Sr[ip+k];
                                    bi = Si[im-k] - Si[ip+k];
                                    cr = Sr[ip-k] - Sr[im+k];
                                    ci = Si[ip-k] - Si[im+k];
                                    tempr+= ar*(br+cr)-ai*(bi-ci);
                                    tempi+= ar*(bi+ci)+ai*(br-cr);
                                }
                            }
                        }else{
                            if (w_flag[mod+u]) {
                                ar = wr[mod+u];
                                ai = wi[mod+u];
                                br = Sr[im];
                                bi = Si[im];
                                tempr+= ar*(br)-ai*(bi);
                                tempi+= ar*(bi)+ai*(br);
                            }
                            for(k=1;k<L+1;k++){
                                if (w_flag[mod+u+k]) {
                                    ar = wr[mod+u+k];
                                    ai = wi[mod+u+k];
                                    br = Sr[im-k];
                                    bi = Si[im-k];
                                    cr = - Sr[im+k];
                                    ci = - Si[im+k];
                                    tempr+= ar*(br+cr)-ai*(bi-ci);
                                    tempi+= ar*(bi+ci)+ai*(br-cr);
                                }
                            }
                        }
                    }
                    r=2;
                    int u = r*(L+1);
                    int im = (m-r)*Np + n;
                    int ip = (m+r)*Np + n;
                    if(r<rframe){
                        if (w_flag[mod+u]) {
                            ar = wr[mod+u];
                            ai = wi[mod+u];
                            br = Sr[im];
                            bi = Si[im];
                            cr = Sr[ip];
                            ci = Si[ip];
                            tempr+= ar*(br+cr)-ai*(bi-ci);
                            tempi+= ar*(bi+ci)+ai*(br-cr);
                        }
                        for(k=1;k<L+1;k++){
                            if (w_flag[mod+u+k]) {
                                ar = wr[mod+u+k];
                                ai = wi[mod+u+k];
                                br = Sr[im-k] + Sr[ip+k];
                                bi = Si[im-k] + Si[ip+k];
                                cr = Sr[ip-k] + Sr[im+k];
                                ci = Si[ip-k] + Si[im+k];
                                tempr+= ar*(br+cr)-ai*(bi-ci);
                                tempi+= ar*(bi+ci)+ai*(br-cr);
                            }
                        }
                    }else{
                        if (w_flag[mod+u]) {
                            ar = wr[mod+u];
                            ai = wi[mod+u];
                            br = Sr[im];
                            bi = Si[im];
                            tempr+= ar*(br)-ai*(bi);
                            tempi+= ar*(bi)+ai*(br);
                        }
                        for(k=1;k<L+1;k++){
                            if (w_flag[mod+u+k]) {
                                ar = wr[mod+u+k];
                                ai = wi[mod+u+k];
                                br = Sr[im-k];
                                bi = Si[im-k];
                                cr = Sr[im+k];
                                ci = Si[im+k];
                                tempr+= ar*(br+cr)-ai*(bi-ci);
                                tempi+= ar*(bi+ci)+ai*(br-cr);
                            }
                        }
                    }
                    
                }else{
                    for(r=1;r<rframe;r++){
                        int u = r*(L+1);
                        int im = (m-r)*Np + n;
                        int ip = (m+r)*Np + n;
                        if (w_flag[mod+u]) {
                            ar = wr[mod+u];
                            ai = wi[mod+u];
                            br = Sr[im];
                            bi = Si[im];
                            cr = Sr[ip];
                            ci = Si[ip];
                            tempr+= ar*(br+cr)-ai*(bi-ci);
                            tempi+= ar*(bi+ci)+ai*(br-cr);
                        }
                        for(k=1;k<L+1;k++){
                            if (w_flag[mod+u+k]) {
                            ar = wr[mod+u+k];
                            ai = wi[mod+u+k];
                            br = Sr[im-k] + Sr[ip+k];
                            bi = Si[im-k] + Si[ip+k];
                            cr = Sr[ip-k] + Sr[im+k];
                            ci = Si[ip-k] + Si[im+k];
                            tempr+= ar*(br+cr)-ai*(bi-ci);
                            tempi+= ar*(bi+ci)+ai*(br-cr);
                            }
                        }
                    }
                    
                    for(r=rframe;r<Q;r++){
                        int u = r*(L+1);
                        int im = (m-r)*Np + n;
                        if (w_flag[mod+u]) {
                            ar = wr[mod+u];
                            ai = wi[mod+u];
                            br = Sr[im];
                            bi = Si[im];
                            tempr+= ar*(br)-ai*(bi);
                            tempi+= ar*(bi)+ai*(br);
                        }
                        for(k=1;k<L+1;k++){
                            if (w_flag[mod+u+k]) {
                                ar = wr[mod+u+k];
                                ai = wi[mod+u+k];
                                br = Sr[im-k];
                                bi = Si[im-k];
                                cr = Sr[im+k];
                                ci = Si[im+k];
                                tempr+= ar*(br+cr)-ai*(bi-ci);
                                tempi+= ar*(bi+ci)+ai*(br-cr);
                            }
                        }
                    }
                }
                
                abstemp = sqrt(pow(tempr, 2.)+pow(tempi, 2.));
                if((abstemp>0)){
                    // update the phase
                    Sr[m*Np+n]= tempr*absspec/abstemp;
                    Si[m*Np+n]= tempi*absspec/abstemp;
                    // propagate changes in the extended regions
                    if((n>=L+1)&&(n<2*L+1)){
                        Sr[m*Np+2*L-n]= Sr[m*Np+n];
                        Si[m*Np+2*L-n]= -Si[m*Np+n];
                    }else if((n>=Nreal-1)&&(n<Naux)){
                        Sr[m*Np+2*Naux-n] = Sr[m*Np+n];
                        Si[m*Np+2*Naux-n] = -Si[m*Np+n];
                    }
                }
            }
        }
    }
}

// General case, slightly slower for Q=2 or Q=4
void Asym_UpdatePhaseanyQ(double *Sr, double *Si, double *wr, double *wi, int *w_flag, double *AmpSpec, int Nreal, int M, int M0, int L, int Q, double threshold, int update)
{
    int n,m,k,r,mod, modneg;
    int Np=Nreal+2*L;
    double tempr, tempi, abstemp, absspec;
    double ar,ai,br,bi,cr,ci;

    int Naux=Nreal+L-1;
    int rframe=0; // indicates the last frame to use on the right
    int cframe=0;

    
    for(m=Q-1;m<(M+Q-1);m++){
        
        rframe= M0+Q-m-1;
        if(rframe>Q){
            rframe=Q;
        }
        cframe=1;
        if(rframe<1){
            cframe=0;
            rframe=1;
        }
        
        for(n=L;n<Naux+1;n++){
            // only update bins whose magnitude is above threshold
            absspec = AmpSpec[m*Np+n];
            if(absspec>threshold){
                tempr=0.;
                tempi=0.;
                mod = ((n-L)%Q)*Q*(L+1);
                modneg = ((Q-((n-L)%Q))%Q)*Q*(L+1);
                if(cframe){
                    if(update==1){
                        tempr+=Sr[m*Np+n]/Q;
                        tempi+=Si[m*Np+n]/Q;
                    }
                    for(k=1;k<L+1;k++){
                        if (w_flag[mod+k]) {
                            ar = wr[mod+k];
                            ai = wi[mod+k];
                            br = Sr[m*Np+n-k];
                            bi = Si[m*Np+n-k];
                            cr = Sr[m*Np+n+k];
                            ci = Si[m*Np+n+k];
                            tempr+= ar*(br+cr)-ai*(bi-ci);
                            tempi+= ar*(bi+ci)+ai*(br-cr);
                        }
                    }
                }

                // contribution of frames m+/-r
                
                for(r=1;r<rframe;r++){ // add frames  at both m+r and m-r (i.e., on left and right)
                    int u = r*(L+1);
                    int im = (m-r)*Np + n;
                    int ip = (m+r)*Np + n;
                    // for frequency n
                    if (w_flag[mod+u]) {
                        ar = wr[mod+u];
                        ai = wi[mod+u];
                        br = Sr[im];
                        bi = Si[im];
                        cr = Sr[ip];
                        ci = Si[ip];
                        tempr+= ar*(br+cr)-ai*(bi-ci);
                        tempi+= ar*(bi+ci)+ai*(br-cr);
                    }
                    // for frequencies n+/-k
                    for(k=1;k<L+1;k++){
                        if (w_flag[mod+u+k]) {
                            ar = wr[mod+u+k];
                            ai = wi[mod+u+k];
                            br = Sr[im-k];
                            bi = Si[im-k];
                            cr = Sr[ip-k];
                            ci = Si[ip-k];
                            tempr+= ar*(br+cr)-ai*(bi-ci);
                            tempi+= ar*(bi+ci)+ai*(br-cr);
                        }
                        if (w_flag[modneg+u+k]) {
                            ar = wr[modneg+u+k];
                            ai = wi[modneg+u+k];
                            br = Sr[ip+k];
                            bi = Si[ip+k];
                            cr = Sr[im+k];
                            ci = Si[im+k];
                            tempr+= ar*(br+cr)-ai*(bi-ci);
                            tempi+= ar*(bi+ci)+ai*(br-cr);
                        }
                    }
                }
                            
                for(r=rframe;r<Q;r++){  // for those only the left frames should be added
                    int u = r*(L+1);
                    int im = (m-r)*Np + n;
                    // for frequency n
                    if (w_flag[mod+u]) {
                        ar = wr[mod+u];
                        ai = wi[mod+u];
                        br = Sr[im];
                        bi = Si[im];
                        tempr+= ar*(br)-ai*(bi);
                        tempi+= ar*(bi)+ai*(br);
                    }
                    // for frequencies n+/-k
                    for(k=1;k<L+1;k++){
                        if (w_flag[mod+u+k]) {
                            ar = wr[mod+u+k];
                            ai = wi[mod+u+k];
                            br = Sr[im-k];
                            bi = Si[im-k];
                            tempr+= ar*(br)-ai*(bi);
                            tempi+= ar*(bi)+ai*(br);
                        }
                        if (w_flag[modneg+u+k]) {
                            ar = wr[modneg+u+k];
                            ai = wi[modneg+u+k];
                            cr = Sr[im+k];
                            ci = Si[im+k];
                            tempr+= ar*(cr)+ai*(ci);
                            tempi+= ar*(ci)-ai*(cr);
                        }
                    }
                }
                
                
                abstemp = sqrt(pow(tempr, 2.)+pow(tempi, 2.));
                if((abstemp>0)){
                    // update the phase
                    Sr[m*Np+n]= tempr*absspec/abstemp;
                    Si[m*Np+n]= tempi*absspec/abstemp;
                    // propagate changes in the extended regions
                    if((n>=L+1)&&(n<2*L+1)){
                        Sr[m*Np+2*L-n]= Sr[m*Np+n];
                        Si[m*Np+2*L-n]= -Si[m*Np+n];
                    }else if((n>=Nreal-1)&&(n<Naux)){
                        Sr[m*Np+2*Naux-n] = Sr[m*Np+n];
                        Si[m*Np+2*Naux-n] = -Si[m*Np+n];
                    }
                }
            }
        }
    }
}


void TF_RTISI_LA(double *Sr, double *Si, double *wr, double *wi, 
        double *wr_asym_init, double *wi_asym_init, double *wr_asym_full, double *wi_asym_full, int *w_flag, int *w_flag_ai, int *w_flag_af,  
        double *AmpSpec, int iter, int LA, int Nreal, int M, int L, int Q, double *ThresholdArray, int update)
{
    int Np=Nreal+2*L;
    double threshold;
    int lframe,nframe;

    for(int m=0;m<M;m++){

        lframe=m-LA;
        nframe=LA;
        if(lframe<0){
            lframe=0;
            nframe=m;
        }
        
        if (Q==2){
            // Initial phase estimate for the newest uncommitted frame
            Asym_UpdatePhaseQ2(Sr+m*Np,Si+m*Np,wr_asym_init,wi_asym_init,w_flag_ai,AmpSpec+m*Np,Nreal,1,0,L,0.,update);
            for(int h=0;h<iter;h++){            
                threshold = ThresholdArray[h];
                // Deal with the other frames using the usual window
                if(LA>0){
                    Asym_UpdatePhaseQ2(Sr+lframe*Np, Si+lframe*Np, wr, wi, w_flag, AmpSpec+lframe*Np, Nreal, nframe, nframe+1, L, threshold, update);
                }
                // Deal with the newest frame first with the asymmetric window
                Asym_UpdatePhaseQ2(Sr+m*Np, Si+m*Np, wr_asym_full, wi_asym_full, w_flag_af, AmpSpec+m*Np, Nreal, 1, 1, L, threshold, update);
            } // end LOOP on iterations
        }else if (Q==4) {
            // Initial phase estimate for the newest uncommitted frame
            Asym_UpdatePhaseQ4(Sr+m*Np, Si+m*Np, wr_asym_init, wi_asym_init, w_flag_ai, AmpSpec+m*Np, Nreal, 1, 0, L, 0., update);
            for(int h=0;h<iter;h++){
                threshold = ThresholdArray[h];
                // Deal with the other frames using the usual window
                if(LA>0){
                    Asym_UpdatePhaseQ4(Sr+lframe*Np, Si+lframe*Np, wr, wi, w_flag, AmpSpec+lframe*Np, Nreal, nframe, nframe+1, L, threshold, update);
                }
                // Deal with the newest frame first with the asymmetric window
                Asym_UpdatePhaseQ4(Sr+m*Np, Si+m*Np, wr_asym_full, wi_asym_full, w_flag_af, AmpSpec+m*Np, Nreal, 1, 1, L, threshold, update);
            } // end LOOP on iterations
        }else {
            // Initial phase estimate for the newest uncommitted frame
        Asym_UpdatePhaseanyQ(Sr+m*Np, Si+m*Np, wr_asym_init, wi_asym_init, w_flag_ai, AmpSpec+m*Np, Nreal, 1, 0, L, Q, 0., update);
            for(int h=0;h<iter;h++){
                threshold = ThresholdArray[h];
                // Deal with the other frames using the usual window
                if(LA>0){
                    Asym_UpdatePhaseanyQ(Sr+lframe*Np, Si+lframe*Np, wr, wi, w_flag, AmpSpec+lframe*Np, Nreal, nframe, nframe+1, L, Q, threshold, update);
                }
                // Deal with the newest frame first with the asymmetric window
                Asym_UpdatePhaseanyQ(Sr+m*Np, Si+m*Np, wr_asym_full, wi_asym_full, w_flag_af, AmpSpec+m*Np, Nreal, 1, 1, L, Q, threshold, update);
            } // end LOOP on iterations

        }
    } // end LOOP on frames
}

