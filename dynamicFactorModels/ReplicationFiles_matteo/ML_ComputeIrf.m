% ML_ComputeIrf - Given a set of admissible rotation computes structural IRF & shocks
% 
% [C, sh] = ML_ComputeIrf(B,GR,rsh)
% 
% Inputs:
%   B - Impulse Responses
%  GR - Good Rotations
% rsh - common shocks
% 
% Outputs:
%   C - Structural IRF
%  sh - Structural Shocks
% 

% Written by Mario Forni
% Modified and commented by Matteo Luciani (matteo.luciani@ulb.ac.be)

function [C, sh] = ML_ComputeIrf(B,GR,rsh)

[N q h]=size(B);
gdraws=size(GR,3);                                                          % number of good rotations
   
if gdraws==0; C = ones(N,q,h,0); sh = [];                                   % when no good draws 
else
    C=zeros(N,q,h,gdraws);                                                  % Preallocates
    if nargin == 3; sh=zeros(length(rsh),q,gdraws); end
    for s = 1:gdraws; 
        for k = 1:h; C(:,:,k,s) = B(:,:,k)*GR(:,:,s); end;                  % Structural Impulse Responses
        if nargin == 3, sh(:,:,s) = rsh*GR(:,:,s); else sh = []; end        % Structural Shocks    
    end        
end
   