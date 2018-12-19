% ML_DfmRawImp - Not Identified IRF for Factor Model
% ML_DfmRawImp(x, q, r, p, s,tb)
%   x - Data
%   q - n dynfactors
%   r - n static factors
%   p - n lags VAR Static Factors
%   s - n lags in the MA representation
%  tb - break point(optonal)
% 
% tb is 2x1 first indicate where is the break, second is which subsample
% 

% Written by Mario Forni
% Modified and commented by Matteo Luciani (matteo.luciani@ulb.ac.be)

function [CL, Chi, eta]  = ML_DfmRawImp(x, q, r, p, s,tb)

[T N] = size(x);
WW = diag(std(x));
y = ML_center(x)*(WW^-1);
[F,lambda,chi]=ML_efactors2(y,r,2);                                         % Estimating Static Factors with the normalization lambda'*lamnda/N=I

if (nargin >5)&&(~ isempty(tb));
    [~,~,~, Lambda] = ML_factor_loadings_break(y,F,tb(1));                  % Estimation of factor loadings when there is a breakpoint
    lambda=Lambda(:,:,tb(2));
end

CL=zeros(N,q,s);                                                            % preallocate
Chi = chi*WW + ones(T,1)*mean(x);                                           % Common Component
[A u]=ML_VAR(F,p,0);                                                        % Estimating the Law of Motion fo the Static Factors   
[eta G]=ML_edynfactors2(u,q);                                               % Estimating Dynamic Factors Innovation
B = ML_MA(s,A,0);                                                           % MA Representation of the static factors
for ii=1:s; CL(:,:,ii)=(lambda*B(:,:,ii)*G).*repmat(diag(WW),1,q); end;     % MA Representation of common components

