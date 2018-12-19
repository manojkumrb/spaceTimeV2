% BLL_SLRImpVAR - Computes IRF for I(1) DFM with SLR restrictions (VAR)
% 
% [Imp,sh] = BLL_SLRImpVAR(Y,cd,q,r,p,h,ID,K0,K1)
% 
% Outputs:
%    Imp - Impulse responses
%     Sh - shocks
% 
% Inputs:
%      Y - Data
%     cd - code for data transformation
%      q - n° of common shocks
%      r - n° of common factors
%      p - n° of lags in the VAR
%      h - n° of lags in the MA representation
%     ID - variables on which restrictionms are imposed
%     KO - Short run restrictions
%     K1 - Long Run Restrictions
% 

% Written by Matteo Luciani (matteo.luciani@ulb.ac.be)

function [Imp,sh] = BLL_SLRImpVAR(Y,cd,q,r,p,h,ID,K0,K1,lr)
if nargin==9; lr=0; end
if lr==1; h1=40; else h1=h; end

y=ML_diff(Y);                                                               % 1st Differences
[yy, ~, sy]=ML_Standardize(y);                                              % Standardize
X=ML_detrend(Y,1);                                                            % Detrend Data
Z=X./repmat(sy,size(X,1),1);                                                % BLL - Standardization 
[~, N]=size(y);                                                             % Panel Dimensions
[~,lambda]=ML_efactors2(yy,r,2);                                            % estimate loadings
F=Z*lambda/N;                                                               % BLL Factors
[~,u,~,B]=ML_VAR(F,p,1,h);                                                  % Estimating VAR on Static Factors   
[eta, G]=ML_edynfactors2(u,q);                                              % Estimating Dynamic Factors Innovation
for ii=1:h; CL(:,:,ii)=(lambda*B(:,:,ii)*G).*repmat(sy',1,q); end;          % MA Representation of common components
C = CumImp(CL,cd);                                                          % comulate IRFs if it is necessary
H=ML_SLRrestrictions(ID,C,C(:,:,h1),K0,K1);                                 % Rotation Matrix
Imp = ML_ComputeIrf(C,H);                                                   % Structural IRF
sh=eta*H;                                                                   % Structural Shocks


function CC = CumImp(Imp, Transf)
notransf = find(Transf==0);
firstdiff = find(Transf==1);
seconddiff = find(Transf==2);
CC = Imp*0;
CC(notransf,:,:,:) = Imp(notransf,:,:,:);
CC(firstdiff,:,:,:) = cumsum(Imp(firstdiff,:,:,:),3);
CC(seconddiff,:,:,:) = cumsum(cumsum(Imp(seconddiff,:,:,:),3),3);
