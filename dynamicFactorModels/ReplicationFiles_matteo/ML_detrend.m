% ML_detrend - OLS detrending
%
% [y, yhat, Beta]=ML_detrend(Y);
%   Y - variables to be detrended
%   y - detrended variables 
%   yhat - fitted values
%   Beta - estimated coefficients
%

% Written by Matteo Luciani (matteo.luciani@ulb.ac.be)

function  [y, yhat, Beta]=ML_detrend(Y,jj)
if nargin==1; jj=1; end
[T, N] = size(Y);
cons=ones(T,1); trend=(1:1:T)';

if jj==1;
    x=[cons trend];
else
    x=trend;
    Y=Y-repmat(Y(1,:),T,1);   
end

yhat=NaN(T,N); y=NaN(T,N); Beta=NaN(N,2);                                   % preallocates
for ii=1:N;
    beta=inv(x'*x)*x'*Y(:,ii);                                              % ols estimation
    yhat(:,ii)=x*beta;                                                      % linear trend
    y(:,ii)=Y(:,ii)-yhat(:,ii);                                             % detrended series
    Beta(ii,:)=beta';                                                       % Store estimated coefficients
end;