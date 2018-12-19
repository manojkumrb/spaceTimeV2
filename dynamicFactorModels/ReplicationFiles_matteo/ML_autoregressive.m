%% ML_autoregressive - Estimate an AR(p) model with OLS selecting the lag length with BIC
% 
% [beta resid r2]=ML_autoregressive(y,det,k,p)
%   Input:
%       y = data Matrix (must be a vector)
%       det = 0 noconstant
%       det = 1 constant
%       det = 2 time trend
%       det = 3 constant + time trend
%       k = maximum number of lags
%       p = number of lags (optional)
%
%   Ouputs:
%       beta  = estimated AR coefficient
%       resid = estimated residual
%

% Written by Matteo Luciani (matteo.luciani@ulb.ac.be)

function  [beta, resid, r2, AL, yhat]=ML_autoregressive(y,det,k,p)

t = size(y,1);


%% Selecting the lag lenght with BIC if required %%

if nargin == 3 || isempty(p);
    r2=NaN(k,1); bic=NaN(k,1);                                              % preallocates variables
    for ii=1:1:k;
        xx=NaN(t-ii,ii);                                                    % preallocates variables
        yy(:,1)=y(ii+1:t,1);         
        for j=1:ii; xx(:,j)=y(ii+1-j:t-j,1); end;     
        [beta,u,~,~, r2(ii)]=ML_ols(yy,xx,det);       
        bic(ii) = log(u'*u/t) + (size(beta,1)/t)*log(t);
        clear beta u v esu yy xx
    end;    
    p=0; for ii=1:k; if bic(ii)==min(bic); p=ii; end; end;   
end;


%% Estimate the AR(p) model %%

yy(:,1)=y(p+1:t,1); 
xx=NaN(t-p,p); for j=1:p; xx(:,j)=y(p+1-j:t-j,1); end; 
[beta,resid,~,~,r2]=ML_ols(yy,xx,det);                            
if det==0; g=1; elseif det==2;g=3; else g=2; end;
yhat=yy-resid;
for ii=g:p+g-1; AA(1,1,ii-g+1)=beta(ii,1); end; AL=1; AL(:,:,2:p+1)=-AA;
    


