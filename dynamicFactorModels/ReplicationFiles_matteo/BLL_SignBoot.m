% BLL_SignBoot - Bootstrap Confidence band for I(1) SDFM wih sign restrictions
% 
% [CC,u] = BLL_SignBoot(Y,cd,q,d,r,p,h,nrest,rest,ndraws,horiz,nboot,maxrot,FP)
% 
% Outputs:
%     CC - Impulse responses
%      u - number of boostrap which finds at least one rotation
% 
% Inputs:
%      Y - I(1) Data                                                        (1)
%     cd - code for data transformation                                     (2)
%      q - n° of common shocks
%      d - n° of transitory shocks
%      r - n° of common factors
%      p - n° of lags in the VAR                                            (5)
%      h - n° of lags in the MA representation                              (6)
%  nrest - normalization restrictions                                       (7)
%   rest - restrictions                                                     (8)
% ndraws - n° of draws for the bootstrap algorithm                          (9)
%  horiz - n° of horizon at which the restriction is imposed                (10)
% nrepli - n° of replications                                               (11)
% maxrot - Maximum Number of Rotation                                       (12) optional
%
% Written by Matteo Luciani (matteo.luciani@frb.gov)
% 
% Reference: Barigozzi, M, Lippi, M. and Luciani M. (2016) "Non-Stationary 
%            Dynamic Factor Models for Large Datasets", FEDS 2016-024
% 

function [CC,u, BB] = BLL_SignBoot(Y,cd,q,d,r,p,h,nrest,rest,ndraws,horiz,nboot,maxrot,FP)

if nargin<14; FP=0; end;
if nargin<13; maxrot=ndraws; end;

y=ML_diff(Y);                                                               % 1st Differences
[yy, ~, sy]=ML_Standardize(y);                                              % Standardize
X=ML_detrend(Y);                                                            % Detrend Data
Z=X./repmat(sy,size(X,1),1);                                                % BLL - Standardization 
[T, N]=size(y);                                                             % Panel Dimensions
[~,lambda]=ML_efactors2(yy,r,2);                                            % estimate factors in first differences
F=Z*lambda/N;                                                               % BLL Factors
[AA, ~, ~, ~, ~, v]=BLL_VECM(F,p,r-q+d,h,1);                                % Estimate the VECM
[V, D] = eigs(cov(v),q,'lm'); G=V*(D^.5); u=v*V*(D^-.5);                    % common shocks
vv=u*G';
n=size(u,1);
BaB=0;

%%% Bootstrap After Botstrap Correction %%%
if BaB==1;
    psi=nan(size(AA,1),size(AA,2),nboot);
    fprintf('Computing Correction for Bootsrap After Bootsrap')
    for bb=1:nboot;
        bootsam=ceil(n*rand(n,1));                                              % reshuffling residuals
        Fstar=ML_decomposition(F,AA,vv(bootsam,:),det);                         % Generates Static Factors        
        AAstar=BLL_VECM(Fstar,p,r-q+d,h,1);                                     % Estimate the VECM
        psi(:,:,bb)=AAstar-AA;
    end
    AA2=AA-mean(psi,3);
    disp(' - done')    
    Atilda=ML_VAR_companion_matrix(ML_VAR_polynomial(AA2,det));                 % VAR in companion form
    mu=max(eig(Atilda));                                                        % Max Root of the VAR
    for pp=1:p; AA2(det+1+pp-1:p:end,:)=AA2(det+1+pp-1:p:end,:)/(mu^pp); end    % Correct coefficients so that max root is 1        
else
    AA2=AA;
end



l=0; ll=0;
if FP==1; CC=nan(N,q,h,nboot); else CC=nan(N,q,h,maxrot*nboot); end; BB=CC;             % Preallocates Variables

%%% Bootsrap Algorithm %%%
tic
for bb=1:nboot;        
    bootsam=ceil(n*rand(n,1));                                                          % reshuffling residuals    
    Fstar=ML_decomposition(F,AA2,vv(bootsam,:),1);                                    % Generates Static Factors    
    [~, ~, ~, CLstar, ~, vstar]=BLL_VECM(Fstar,p,r-q+d,h,1);                            % Estimate the VECM
    [V, D] = eigs(cov(vstar),q,'lm'); Gstar=V*(D^.5);                                   % common shocks    
    for jj=1:h; Phistar(:,:,jj)=(lambda*CLstar(:,:,jj)*Gstar).*repmat(sy',1,q); end;	% Impulse Responses
    PHIstar = CumImp(Phistar,cd);                                                       % comulate IRFs if it is necessary    
    GR = ML_GoodPolar(PHIstar,nrest,rest,ndraws,horiz,maxrot);                          % Good Rotations
    C = ML_ComputeIrf(PHIstar,GR); B = ML_ComputeIrf(Phistar,GR);                       % Compute IRF
            
    ngr=size(GR,3);
    if ngr ~= 0; 
        l = ll+1; ll=ll+ngr;                                                            % count how many bootstrap yeld at least one rotation              
        if (ngr > 1) && (FP==1);
            BR=ML_SignRestrictionFryPagan(C,rest);                                      % Identify the best rotation in the sense of Fry & Pagan                                           
            C=C(:,:,:,BR); B=B(:,:,:,BR); ll=l;                                         % Keep Only the FP rotation     
        end  
        CC(:,:,:,l:ll) = C; BB(:,:,:,l:ll) = B;                                         % Store IRF     
    end;         
    disp([nboot-bb ngr toc])                                                            % Counter
    tic
end

CC(:,:,:,ll+1:end)=[];
BB(:,:,:,ll+1:end)=[];
disp([ll l]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CC = CumImp(Imp, Transf)
notransf = find(Transf==0);
firstdiff = find(Transf==1);
seconddiff = find(Transf==2);
CC = Imp*0;
CC(notransf,:,:,:) = Imp(notransf,:,:,:);
CC(firstdiff,:,:,:) = cumsum(Imp(firstdiff,:,:,:),3);
CC(seconddiff,:,:,:) = cumsum(cumsum(Imp(seconddiff,:,:,:),3),3);




