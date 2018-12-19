% ML_VAR - Estimates a VAR(k) for a vector of variables y
% if y is a single variable it estimates an AR model with OLS
%
% [A,u,AL,CL]=ML_VAR(y,k,jj);
%   y  = vector of endogenous variables
%   k  = number of lags
%   jj = 0 noconstant
%   jj = 1 constant
%   jj = 2 time trend
%   jj = 3 constant + time trend
%   A  = Matrix of coefficients for the reduced form 
%   u  = Vector of Residuals
%
%  y_t=A(L)y_{t-1}+u_t
%  y_t=C(L)u_t
% ----------------------------------------------------------

% written by Matteo Luciani (matteo.luciani@ulb.ac.be)

function [A,u,AL,CL,C1]=ML_VAR(y,k,jj,s)
[T N] = size(y);

%%% Building Up the vector for OLS regression %%%
for i=1:N,
    yy(:,i)=y(k+1:T,i);
    for j=1:k,
        xx(:,k*(i-1)+j)=y(k+1-j:T-j,i);
    end;
end;


%%% OLS Equation-By-Equation %%%
for i=1:N;
    [A(:,i),u(:,i)]=ML_ols(yy(:,i),xx,jj);
end;


At=A; if jj==3; At(1:2,:)=[]; elseif jj==1||jj==2; At(1,:)=[]; end
AL=NaN(N,N,k); for kk=1:k; AL(:,:,kk)=At(1+kk-1:k:end,:)'; end  


%%% Impulse Responses %%%
if nargin<4;s=20; end
CL(:,:,1) = eye(N);
for ss=2:s;
    CL(:,:,ss) = 0;
    for ii = 1:min(ss-1,k);        
        temp3=AL(:,:,ii)*CL(:,:,ss-ii);        
        CL(:,:,ss)=CL(:,:,ss)+temp3;
    end;
end;

C1=eye(N); for i=1:k; C1=C1-AL(:,:,i); end; C1=inv(C1);                     % Long Run Multipliers


