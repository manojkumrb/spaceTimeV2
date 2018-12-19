% BLL_SignVAR - IRF sign restrictions VAR Level estimated with BLL
%               same as BLL_Sign but it estmates a VAR in Levels for Ft
% 
% [C,e,B]=BLL_SignVAR(Y,cd,q,r,p,h,nrest,rest,ndraws,horiz,maxrot,FP)
% 
% Outputs:
%      C - Impulse responses
%      e - residuals
% 
% Inputs:
%      Y - Data                                                             (1)
%     cd - code for data transformation                                     (2)
%      q - n° of dynamic factors                                            (3)
%      r - n° of static factors                                             (4)
%      p - n° of lags in the VAR                                            (5)
%      h - n° of lags in the MA representation                              (6)
%  nrest - normalization restrictions                                       (7)
%   rest - restrictions                                                     (8)
% ndraws - n° of draws for the bootstrap algorithm                          (9)
%  horiz - n° of horizon at which the restriction is imposed                (10)
% nrepli - n° of replications                                               (11)
% maxrot - Maximum Number of Rotation                                       (12) optional
%     FP - Fry and Pagan Correction
%
% Written by Matteo Luciani (matteo.luciani@frb.gov)
% 
% Reference: Barigozzi, M, Lippi, M. and Luciani M. (2016) "Non-Stationary 
%            Dynamic Factor Models for Large Datasets", FEDS 2016-024
% 

function [C,e,B] = BLL_SignVAR(Y,cd,q,r,p,h,nrest,rest,ndraws,horiz,maxrot,FP)

if nargin<12; FP=0; end;
if nargin<11; maxrot=ndraws; end;

y=ML_diff(Y);                                                               % 1st Differences
[yy, ~, sy]=ML_Standardize(y);                                              % Standardize
X=ML_detrend(Y);                                                            % Detrend Data
Z=X./repmat(sy,size(X,1),1);                                                % BLL - Standardization 

[~, N]=size(y);                                                             % Panel Dimensions
[~,lambda]=ML_efactors2(yy,r,2);                                            % estimate factors in first differences
F=Z*lambda/N;                                                               % BLL Factors
[~,v,~,CL]=ML_VAR(F,p,1,h);                                                 % Estimating VAR on Factors   
[V D] = eigs(cov(v),q,'lm'); G=V*(D^.5); u=v*V*(D^-.5);                     % common shocks
for jj=1:h; Phi(:,:,jj)=(lambda*CL(:,:,jj)*G).*repmat(sy',1,q); end;        % Impulse Responses
CC = CumImp(Phi,cd);                                                        % comulate IRFs if it is necessary

GR = ML_GoodPolar(CC,nrest,rest,ndraws,horiz,maxrot);                       % Good Rotations
[C e] = ML_ComputeIrf(CC,GR,u); B = ML_ComputeIrf(Phi,GR);                  % Compute IRF

if size(C,4)==0;    
else    
    if FP==1;                                                               % In case the Fry Pagan Correction is Required
        BR=ML_SignRestrictionFryPagan(C,rest);                              % Identify the best rotation in the sense of Fry & Pagan
        C=C(:,:,:,BR); e=e(:,:,BR); B=B(:,:,:,BR);                          % Keep Only the FP rotation
    end
end



% 
% [CC,Y,u] = ML_SignBoot(X,cd,q,r,p,h,L,nrest,rest,ndraws,horiz,maxrot)
% 
% Outputs:
%     CC - Impulse responses
%      Y - Percentiles of the distribution of the IRF 
%      u - number of boostrap which finds at least one rotation
% 
% Inputs:
%      X - Data                                                             (1)
%     cd - code for data transformation                                     (2)
%      q - n° of dynamic factors                                            (3)
%      r - n° of stati factors                                              (4)
%      p - n° of lags in the VAR                                            (5)
%      h - n° of lags in the MA representation                              (6)
%  nrest - normalization restrictions                                       (7)
%   rest - restrictions                                                     (8)
% ndraws - n° of draws for the bootstrap algorithm                          (9)
%  horiz - n° of horizon at which the restriction is imposed                (10)
% nrepli - n° of replications                                               (11)
% maxrot - Maximum Number of Rotation                                       (12) optional
%
% Written by Mario Forni. 
% Modified and Commented by Matteo Luciani (matteo.luciani@ulb.ac.be)
% 
% Reference: Matteo Luciani (2012) "Monetary Policy and the Housing Market: 
%            A Structural Factor Analysis", Ecares Working Paper 2012-035
%            to appear on the Journal of Applied Econometrics.
% 




