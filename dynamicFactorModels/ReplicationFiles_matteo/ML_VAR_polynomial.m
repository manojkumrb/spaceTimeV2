% ML_VAR_polynomial - Construct the Autoregressive Polynomial (RHS)
%
% Y(t)=AL(:,:,1)*Y(t-1)+...+AL(:,:,p)*Y(t-p)+u(t)
%
% AL=ML_VAR_polynomial(A,det)
%   A   - OLS Estimates Obtained with ML_VAR
%   det - Deterministic Component
%

% Written by Matteo Luciani (matteo.luciani@ulb.ac.be)

function AL=ML_VAR_polynomial(A,det)

[g N]=size(A);
%%% Retrieving the Number of Lags in the VAR %%%
if      det==0; p=(g/N);            
elseif  det==3; p=((g-2)/N); 
else            p=((g-1)/N);  end

if      det==0; temp1 = A; 
elseif  det==3; temp1 = A(3:g,:);
else            temp1 = A(2:g,:); end
temp2=temp1'; 
AL= zeros(N, N, p);
for i=1:p; m=i; for j=1:N; AL(:,j,i)= temp2(:,m); m=m+p;  end; end;