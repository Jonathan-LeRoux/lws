function spec = stft(x,flength,fshift,w,opts)

% Short-Time Fourier Transform for single or multichannel signals
% Jonathan Le Roux.

% Inputs:
%  x            data vector (nsamples x nchannels)
%  flength      frame length 
%  fshift       frame shift 
%  w            analysis window function (default: sqrt(hanning))
%  opts         optional arguments
%               opts.framepadding: if 1, 0-pad signal before applying STFT
%
% Output:
%  spec         STFT spectrogram (nbin x nframes x nchannels)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2008-2017 Jonathan Le Roux
%   Apache 2.0  (http://www.apache.org/licenses/LICENSE-2.0) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if mod(flength,2) == 1
    error('odd ffts not implemented')
end

if size(x,2)>size(x,1)
    x=x'; % ensure the channels correspond to columns of x
end

if nargin < 3
   fshift=flength/2; 
end
if ~exist('opts','var')
    opts = struct;
end
if ~isfield(opts, 'framepadding')
    opts.framepadding = 0;
end 

I= size(x,2);
Q=flength/fshift;
if opts.framepadding == 1 
    %padding with 0 for perfect reconstruction near the boundaries
    x=cat(1,zeros((Q-1)*fshift,I),x,zeros((Q-1)*fshift,I));
end
T = size(x,1);

M=ceil((T-flength)/fshift)+1;
% Make sure we have an entire frame at the end
x=cat(1,x,zeros((M-1)*fshift+flength-T,I));

spec=zeros(flength,M,I);
if nargin<4
    w=sqrt((0.5-0.5*cos(2*pi*(1:2:2*flength-1)'/(2*flength)))/Q*2);
else
    if size(w,2)>1, w=w'; end
    if length(w)~=flength
        error('The size of the specified window is incorrect'); 
    end
end

for m=1:M
   frame=bsxfun(@times,x((1:flength)+(m-1)*fshift,:),w);
   temp=fft(frame);% we use the same normalization as matlab, i.e. 1/T in ifft only
   spec(:,m,:)=temp;
end

spec=spec(1:(flength/2+1),:,:);