% Example script for phase reconstruction using Local Weighted Sums (LWS)
%
%
% If you use this code, please cite:
% ## LWS ##
% Jonathan Le Roux, Hirokazu Kameoka, Nobutaka Ono and Shigeki Sagayama, 
% "Fast Signal Reconstruction from Magnitude STFT Spectrogram Based on Spectrogram Consistency," 
% in Proc. International Conference on Digital Audio Effects (DAFx), pp. 397--403, Sep. 2010.
% @InProceedings{LeRoux2010DAFx09,
%   author =	 {Jonathan {Le Roux} and Hirokazu Kameoka and Nobutaka Ono and Shigeki Sagayama},
%   title =	 {Fast Signal Reconstruction from Magnitude {STFT} Spectrogram Based on Spectrogram Consistency},
%   booktitle =	 {Proc. International Conference on Digital Audio Effects (DAFx)},
%   year =	 2010,
%   pages =	 {397--403},
%   month =	 sep
% }
% ## Online LWS, "No future" LWS ##
% Jonathan Le Roux, Hirokazu Kameoka, Nobutaka Ono and Shigeki Sagayama, 
% "Phase initialization schemes for faster spectrogram-consistency-based signal reconstruction," 
% Proc. of ASJ Autumn Meeting, 3-10-3, Sep. 2010.
% @InProceedings{LeRoux2010ASJ09,
%   author =	 {Jonathan {Le Roux} and Hirokazu Kameoka and Nobutaka Ono and Shigeki Sagayama},
%   title =	 {Phase Initialization Schemes for Faster Spectrogram-Consistency-Based Signal Reconstruction},
%   year =	 2010,
%   booktitle =	 {Proceedings of the Acoustical Society of Japan Autumn Meeting (ASJ)},
%   number =	 {3-10-3},
%   month =	 mar
% }
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2008-2017 Jonathan Le Roux
%   Apache 2.0  (http://www.apache.org/licenses/LICENSE-2.0) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with the wavread/audioread and wavwrite/audiowrite annoyance
useaudioread = 0;
if exist('audioread','file')
    useaudioread = 1;
end
%% Pipeline
do_no_future_lws = 1;
do_online_lws    = 1;
do_batch_lws     = 1;

%% Wav output flag
output_wav = 1;
%% Set analysis parameters
N=512;
wshift=128;
Q=N/wshift;
% Note: the code assumes the windows are symmetric, i.e., W(n) = W(N+1-n)
%W=sqrt((0.5-0.5*cos(2*pi*(1:2:2*N-1)'/(2*N)))/Q*2);
%S=sqrt((0.5-0.5*cos(2*pi*(1:2:2*N-1)'/(2*N)))/Q*2);
W=sqrt((0.5-0.5*cos(2*pi*(0:(N-1))'/(N)))/Q*2);
S=sqrt((0.5-0.5*cos(2*pi*(0:(N-1))'/(N)))/Q*2);

%% Get the wav file in
infile='test.wav'; % please provide a test file
if useaudioread
    [x,fs]=audioread(infile);
else
    [x,fs]=wavread(infile);
end
if size(x,2)>1, fprintf('This code only handles single channel files\n'); return; end
fprintf(1,'Processing file: %s (length: %.2f s)\n',infile,length(x)/fs);
X=stft(x,N,wshift,W);
Xpow = sum(abs(X(:)).^2);

tmp = stft(istft(abs(X),wshift,S),N,wshift,W)-X;
C = 10*log10(Xpow / sum(abs(tmp(:)).^2));
fprintf(1,'Abs(X)          : %5.2f dB\n',C);

%% Create weights
L=5;
weights=create_weights(W,S,wshift,L);

%% Create weights for asymmetric windows
[W_asym_init,W_asym_full] = build_asymmetric_windows(W.*S,wshift);
weights_asym_init=create_weights(W_asym_init,S,wshift,L);
weights_asym_full=create_weights(W_asym_full,S,wshift,L);

%% Perform "no future" initialization
% At this point, this does not seem to help much w.r.t. doing nothing
if do_no_future_lws
    tic;
    X0=nofuture_lws(X,weights_asym_init,0);
    time_nofuture_lws = toc;
    
    tmp = stft(istft(X0,wshift,S),N,wshift,W)-X0;
    C0 = 10*log10(Xpow / sum(abs(tmp(:)).^2));
    fprintf(1,'No future init  : %5.2f dB (time: %.2f s)\n',C0,time_nofuture_lws);
    if output_wav
        x0=istft(X0,wshift,S);
        if useaudioread
            audiowrite('out_nofuture_lws.wav',x0,fs);
        else
            wavwrite(x0,fs,'out_nofuture_lws.wav');
        end
    end

else
    X0 = abs(X);
end

%% Perform a few iterations of Online LWS (aka. TF-RTISI-LA)
if do_online_lws
    tic;
    look_ahead        = 3;
    online_iterations = 10;
    online_alpha      = 1;
    online_beta       = 0.1;
    online_gamma      = 1;
    online_thresholds = online_alpha*exp(-online_beta*(0:(online_iterations-1)).^online_gamma);
    
    X1=online_lws(X0,weights,weights_asym_init,weights_asym_full,online_thresholds,look_ahead);
    time_online_lws = toc;
    
    tmp = stft(istft(X1,wshift,S),N,wshift,W)-X1;
    C1 = 10*log10(Xpow / sum(abs(tmp(:)).^2));
    fprintf(1,'+Online LWS     : %5.2f dB (time: %.2f s)\n',C1,time_online_lws);
    if output_wav
        x1=istft(X1,wshift,S);
        if useaudioread
            audiowrite('out_online_lws.wav',x1,fs);
        else
            wavwrite(x1,fs,'out_online_lws.wav');
        end
    end

else
    X1 = X0;
end

%% run batch LWS
if do_batch_lws
    tic;
    iterations = 100;
    alpha      = 100;
    beta       = 0.1;
    gamma      = 1;
    thresholds = alpha*exp(-beta*(0:(iterations-1)).^gamma);

    Y=batch_lws(X1,weights,thresholds);
    time_bach_lws = toc;
    
    tmp = stft(istft(Y,wshift,S),N,wshift,W)-Y;
    Cfinal = 10*log10(Xpow / sum(abs(tmp(:)).^2));
    fprintf(1,'+Batch LWS      : %5.2f dB (time: %.2f s)\n',Cfinal,time_bach_lws);
    if output_wav
        y=istft(Y,wshift,S);
        if useaudioread
            audiowrite('out_batch_lws.wav',y,fs);
        else
            wavwrite(y,fs,'out_batch_lws.wav');
        end
    end

else
    Y = X1;
end

