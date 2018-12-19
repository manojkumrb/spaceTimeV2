% ML_DfmSLRImp - Structural Factor Model Estimation with Zero Restrictions
%
% [Imp, Dimp, sh] = ML_DfmSLRImp(x,cd,q,r,p,h,K0,K1)
%
%       sh - shocks
% 
%        x - Data                                                           (1)
%       cd - order of differentiation of the variables                      (2)
%        q - number of dynamic factors                                      (3)
%        r - number of static factors                                       (4)
%        p - number of lags in the VAR for the static factors               (5)
%        h - number of lags of the Impulse Responses                        (6)
%       ID - variables on which restrictions are imposed
%       K0 - shor run restrictions
%       K1 - Long Run restrictions
% 

% Written by Matteo Luciani (matteo.luciani@ulb.ac.be)

function [Imp,Dimp,sh] = ML_DfmSLRImp(x,cd,q,r,p,h,ID,K0,K1)

[~, N] = size(x);
CL=zeros(N,q,h); Dimp=CL; Imp=Dimp;                                         % preallocate
WW = diag(std(x));                                                          % standard deviation
y = ML_center(x)*(WW^-1);                                                   % standardize
[F,lambda]=ML_efactors2(y,r,2);                                             % Estimating Static Factors with the normalization lambda'*lamnda/N=I
[~,u,~,B]=ML_VAR(F,p,0,h);                                                  % Estimating VAR on Static Factors   
[eta G]=ML_edynfactors2(u,q);                                               % Estimating Dynamic Factors Innovation
for ii=1:h; CL(:,:,ii)=(lambda*B(:,:,ii)*G).*repmat(diag(WW),1,q); end;     % MA Representation of common components
C = CumImp(CL,cd);                                                          % comulate IRFs if it is necessary
H=ML_SLRrestrictions(ID,C,C(:,:,end),K0,K1);                                % Rotation Matrix
Dimp = ML_ComputeIrf(CL,H); Imp = ML_ComputeIrf(C,H);                       % Structural IRF
sh=eta*H;                                                                   % Structural Shocks


