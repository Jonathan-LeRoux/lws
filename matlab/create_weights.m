function weights = create_weights(W,S,fshift,L)

% Compute the (L+1)xQxQ matrix of complex weights used by the LWS code
%
% Inputs:
%  w            analysis window function
%  s            synthesis window function
%  fshift       frame shift
%  L            truncation order (in frequency direction) in LWS
%
% Output:
%  weights      (L+1)xQxQ matrix of complex weights
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2008-2017 Jonathan Le Roux
%   Apache 2.0  (http://www.apache.org/licenses/LICENSE-2.0) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(W,1)>1
    W=W';
end
if size(S,1)>1
    S=S';
end

T=length(W);
interval=0:(L);
expinterv=exp(-1i*2*pi*interval'*(0:(T-1))/T);
Q=T/fshift;
windowprod=zeros(T,Q);
for q=1:Q
    index=1:(T-(q-1)*fshift);
    windowprod(index,q)=W(index).*S(index+(q-1)*fshift)/T;
end
weights=(expinterv*windowprod).*exp(-1i*2*pi*interval'*(0:(Q-1))/Q);
weights(1,1) = weights(1,1) - 1;

tmp=shiftdim(exp(1j*2*pi*(0:(Q-1))'*(0:(Q-1))/Q),-1);
weights=bsxfun(@times,weights,tmp);