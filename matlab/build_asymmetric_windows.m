function [W_asym_init,W_asym_full] = build_asymmetric_windows(WS,fshift)

% This computes the mirrored envelope used in Zhu et al.'s RTISI-LA
% Note that the input WS should be the *product* of the analysis and 
% synthesis windows.
%
% Inputs:
%  WS           analysis window function
%  fshift       frame shift
%
% Output:
%  W_asym_init  assymetric window used for initialization in RTISI-LA
%  W_asym_full  assymetric window used for the newest frame in RTISI-LA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2008-2017 Jonathan Le Roux
%   Apache 2.0  (http://www.apache.org/licenses/LICENSE-2.0) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=length(WS);
Q=N/fshift;

tmp = zeros(N,Q);
for q=0:(Q-1)
    tmp(1:(N-q*fshift),q+1)=WS((q*fshift+1):N);
end

W_asym_init = flipud(sum(tmp(:,2:Q),2));
W_asym_full = flipud(sum(tmp,2));

% for Q=2, we fall back to W as there is no overlap between the incomplete
% window and its symmetric w.r.t. to the center of the interval
if Q==2
   W_asym_init = WS;    
end

end
