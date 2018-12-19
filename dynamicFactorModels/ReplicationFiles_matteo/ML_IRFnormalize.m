% ML_IRFnormalize - Normalize Impulse Responses
%
% IRF=ML_IRFnormalize(IRF,norm,nV)
% 
%    IRF - s by n by q matrix of Impulse Responses
%   norm - 1 by q vector containing the normalization restrictions
%     nV - location of the variable on which the normalization need to be applied
% 

% written by Matteo Luciani (matteo.luciani@ulb.ac.be)

function [IRF irf]=ML_IRFnormalize(IRF,norm,nV,irf)

if nargin==4; ev=1; else ev=0;irf=[]; end

q=size(IRF,3);

for qq=1:q;
    a=norm(qq)/IRF(1,nV(qq),qq);
    IRF(:,:,qq)=a*IRF(:,:,qq);
    if ev==1; irf(:,:,qq)=a*irf(:,:,qq); end;
end