function signal = istft(spec,fshift,s,opts)

% Inverse STFT by overlap-add for single or multichannel signals
% Jonathan Le Roux
%
% Inputs:
%  spec         STFT spectrogram (nbin x nframes x nchannels)
%  fshift       frame shift
%  s            synthesis window function (default: sqrt(hanning))
%  opts         optional arguments
%               opts.framepadding: indicates if the signal was 0-padded
%                                  before applying STFT
%
% Output:
%  signal       resynthesized signal (each channel is a column vector)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2008-2017 Jonathan Le Roux
%   Apache 2.0  (http://www.apache.org/licenses/LICENSE-2.0) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N,M,I]=size(spec);
if mod(N,2)==1
    flength=2*(N-1);
else
    flength=N;
end
if ~exist('opts','var')
    opts = struct;
end

if ~isfield(opts, 'framepadding')
    opts.framepadding = 0;
end 

Q=flength/fshift;
if ~exist('s','var')||isempty(s)
    s=sqrt((0.5-0.5*cos(2*pi*(1:2:2*flength-1)'/(2*flength)))/Q*2);
else
    if size(s,2)>1, s=s'; end
    if length(s)~=flength
        error('The size of the specified window is incorrect'); 
    end
end

T=fshift*(M-1)+flength;
signal=zeros(T,I);

for m=1:M
   iframe=ifft(reshape(spec(:,m,:),[],I),flength,'symmetric');
   signal((m-1)*fshift+(1:flength),:)=signal((m-1)*fshift+(1:flength),:)...
       + bsxfun(@times,iframe,s);    
end

if opts.framepadding == 1
    signal=signal(((Q-1)*fshift+1):(T-(Q-1)*fshift),:);    
end