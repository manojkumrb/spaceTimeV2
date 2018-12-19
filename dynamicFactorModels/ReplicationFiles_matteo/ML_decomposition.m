% ML_decomposition - Given a series of shocks/residuald u construct Yhat=inv(I-A(L))*u=B(L)u
%
% yhat=ML_decomposition(y,beta,u,det,H)
%   yhd T by q by n
%   beta coefficients of the VAR, p*n+det by n
%   u = structural shocks
%   det = deterministic component in the VAR
%   H rotation matrix for structural shocks (optional)
%   v = u*H
%

% Written by Matteo Luciani (matteo.luciani@ulb.ac.be)

function yhat=ML_decomposition(y,beta,u,det,H)

[g N]=size(beta); T=size(y,1); yhat=zeros(T,N); 
if nargin<5; H=eye(N,N); end;  q=size(u,2);

%%% Retrieving the Number of Lags in the VAR %%%
if      det==0; p=(g/N);            
elseif  det==3; p=((g-2)/N); 
else            p=((g-1)/N); ; end

%%% Setting-up the deterministic component %%%
if      det==1; a=beta(1,:);
elseif  det==2; b=beta(1,:);
elseif  det==3; a=beta(1,:); b=beta(2,:); end

A=ML_VAR_polynomial(beta,det);

%%% Starting Values %%%
for i=1:p; yhat(i,:)=y(i,:); end;

%%% Rolling Over the model %%%
u=[zeros(p,q);u];
for t=p+1:T;    
    for ii=1:p;            
        for jj=1:N;            
            yhat(t,jj)=yhat(t,jj)+A(jj,:,ii)*yhat(t-ii,:)';            
        end        
    end;                
    if     det==0; yhat(t,:)=yhat(t,:)      +(u(t,:)*H);        
    elseif det==1; yhat(t,:)=a+yhat(t,:)    +(u(t,:)*H);
    elseif det==2; yhat(t,:)=b*t+yhat(t,:)  +(u(t,:)*H);
    elseif det==3; yhat(t,:)=a+b*t+yhat(t,:)+(u(t,:)*H);
    end   
end    


