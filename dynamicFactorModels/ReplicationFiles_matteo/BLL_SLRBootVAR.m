% BLL_SignBootVAR - Bootstrap Confidence band for I(1) SDFM wih sign restrictions
% 
% [CC,u] = BLL_SignBoot(Y,cd,q,d,r,p,h,nrest,rest,ndraws,horiz,nboot,maxrot,FP)
% 
% Outputs:
%     CC - Impulse responses
%      u - number of boostrap which finds at least one rotation
% 
% Inputs:
%      Y - I(1) Data                                                             (1)
%     cd - code for data transformation                                     (2)
%      q - n° of common shocks
%      d - n° of transitory shocks
%      r - n° of common factors
%      p - n° of lags in the VAR                                            (5)
%      h - n° of lags in the MA representation                              (6)
%     ID - variables on which restrictionms are imposed
%     KO - Short run restrictions
%     K1 - Long Run Restrictions

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

function [CC,eflag] = BLL_SLRBootVAR(Y,cd,q,r,p,h,ID,K0,K1,nboot,lr)

if nargin==10; lr=0; end
if lr==1; h1=40; else h1=h; end
h1

y=ML_diff(Y);                                                               % 1st Differences
[yy, ~, sy]=ML_Standardize(y);                                              % Standardize
X=ML_detrend(Y,1);                                                            % Detrend Data
Z=X./repmat(sy,size(X,1),1);                                                % BLL - Standardization 
[~, N]=size(y);                                                             % Panel Dimensions
[~,lambda]=ML_efactors2(yy,r,2);                                            % estimate factors in first differences
F=Z*lambda/N;                                                               % BLL Factors
[AA,v,~,B]=ML_VAR(F,p,1,h);                                                 % Estimating VAR on Factors   
[V, D] = eigs(cov(v),q,'lm'); G=V*(D^.5); u=v*V*(D^-.5);                    % common shocks
vv=u*G';                                                                    % rxr residuals with rank q
n=size(u,1);
det=1;

%%% Bootstrap After Botstrap Correction %%%
BaB=0;
if BaB==1;
    psi=nan(size(AA,1),size(AA,2),nboot);
    fprintf('Computing Correction for Bootsrap After Bootsrap')
    for bb=1:nboot;
        bootsam=ceil(n*rand(n,1));                                          % reshuffling residuals
        Fstar=ML_decomposition(F,AA,vv(bootsam,:),det);                     % Generates Static Factors                       
        AAstar=ML_VAR(Fstar,p,1,h);                                         % Estimating VAR on Factors
        
        psi(:,:,bb)=AAstar-AA;
    end
    AA2=AA-mean(psi,3);
    disp(' - done')    
%     Atilda=ML_VAR_companion_matrix(ML_VAR_polynomial(AA2,det),det);             % VAR in companion form
%     mu=max(eig(Atilda));                                                        % Max Root of the VAR
%     for pp=1:p; AA2(det+1+pp-1:p:end,:)=AA2(det+1+pp-1:p:end,:)/(mu^pp); end    % Correct coefficients so that max root is 1        
else
    AA2=AA;
end

%%% Bootsrap Algorithm %%%
tic
fprintf('Running Bootsrap Algorithm')
WW=repmat(sy',1,q);
CC=nan(N,q,h,nboot);                                                        % Preallocates Variables
for bb=1:nboot;        
    bootsam=ceil(n*rand(n,1));                                              % reshuffling residuals    
    Fstar=ML_decomposition(F,AA2,vv(bootsam,:),det);                        % Generates Static Factors    
    [AAstar,vstar,~,CLstar]=ML_VAR(Fstar,p,1,h);                                                 % Estimating VAR on Factors           
    [~, Gstar]=ML_edynfactors2(vstar,q);                                    % common shocks       
    for jj=1:h; Phistar(:,:,jj)=(lambda*CLstar(:,:,jj)*Gstar).*WW; end;     % Impulse Responses
    BB = CumImp(Phistar,cd);                                                % comulate IRFs if it is necessary        
    [H eflag(:,bb)]=ML_SLRrestrictions(ID,BB,BB(:,:,h1),K0,K1);                          % Rotation Matrix    
%     BB2= ML_ComputeIrf(BB,H);         
%     plot([squeeze(BB(31,1,:)) squeeze(BB2(31,1,:))])

    CC(:,:,:,bb) = ML_ComputeIrf(BB,H);                                     % Compute IRF                       
end
disp(' - done')
toc



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




