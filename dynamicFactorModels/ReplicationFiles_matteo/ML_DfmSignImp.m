% ML_DfmSignImp - Structural Factor Model Estimation with Sign Restrictions
%
% [Imp, Dimp, sh] = DfmSignImp(X,cd, q, r, k, h, nrest, rest, ndraws, horiz)
%
%       sh - shocks
% 8
%        X - Data                                                           (1)
%       cd - order of differentiation of the variables                      (2)
%        q - number of dynamic factors                                      (3)
%        r - number of static factors                                       (4)
%        k - number of lags in the VAR for the static factors               (5)
%        h - number of lags of the Impulse Responses                        (7)
%    nrest - normalization restrictions                                     (8)
%     rest - cell containing the sign restrictions                          (9)
%   ndraws - number of draws for the bootstrap algorithm                    (10)
%    horiz - horizon in which sign restrictions are imposed                 (11) optional
%       FP - if 1 uses the Fry & Pagan Modification                         (12) optional
% 

% Written by Mario Forni
% Modified and commented by Matteo Luciani (matteo.luciani@ulb.ac.be)

function [Imp,Dimp,sh,imp] = ML_DfmSignImp(X,cd,q,r,k,h,nrest,rest,ndraws,horiz,maxrot,FP,tb)

if nargin<13; tb=[]; end;
if nargin<12; FP=0; end;
if nargin<11; maxrot=ndraws; end;

[B, ~, rsh]  = ML_DfmRawImp(X, q, r, k, h,tb);                              % Compute impulse response for non structural shocks
C = CumImp(B,cd);                                                           % comulate IRFs whether it is necessary
GR = ML_GoodPolar(C, nrest, rest, ndraws, horiz,maxrot);                    % look for admissible rotations
[Imp, sh] = ML_ComputeIrf(C,GR,rsh);                                        % computes structural IRF & shocks for all admissible rotation
Dimp = ML_ComputeIrf(B,GR);                                                 % computes structural IRF & shocks for all admissible rotation


if size(Imp,4)==0;
imp=[];
else
    imp=Imp;
if FP==1;                                                                   % In case the Fry Pagan Correction is Required
    BR=ML_SignRestrictionFryPagan(Imp,rest);                                % Identify the best rotation in the sense of Fry & Pagan           
    Imp=Imp(:,:,:,BR); Dimp=Dimp(:,:,:,BR); sh = sh(:,:,BR);                % Keep Only the FP rotation
end 
end
