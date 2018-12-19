% ML_SLRBoot - Bootstrap Confidence band for SDFM wih zero restrictions
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

% Written by by Matteo Luciani (matteo.luciani@ulb.ac.be)

function [CC eflag] = ML_SLRBoot(x,cd,q,r,p,h,ID,K0,K1,nboot)

N = size(x,2);
WW = diag(std(x)); WWW=repmat(diag(WW),1,q);                                % Standard deviation
y = ML_center(x)*(WW^-1);                                                   % standardize variables
det=1;                                                                      % by default a constant is included
[fhat,lambda]=ML_efactors2(y,r,2);                                          % Estimating Static Factors with the normalization lambda'*lamnda/N=I
[A u]=ML_VAR(fhat,p,det);                                                   % Estimating the Law of Motion fo the Static Factors   
[eta G]=ML_edynfactors2(u,q);                                               % Estimating Dynamic Factors Innovation
A2=ML_VARbias(fhat,u,A,p,det,nboot,eta,G');                                 % Computes the bias induced by VAR Estimation - Kilian (1998)
% A2=A;
n=size(u,1);


disp('Computing Bootstrap Confidence band - Please Wait')
CC=nan(N,q,h,nboot);                                                        % Preallocates Variables
tic
for bb=1:nboot; 
    bootsam=ceil(n*rand(n,1));                                              % reshuffling residuals    
    fstar=ML_decomposition(fhat,A2,eta(bootsam,:),det,G');                  % Generates Static Factors        
    [~,ustar,~,Bstar]=ML_VAR(fstar,p,det,h);                                % Estimating VAR on Static Factors       
    [~, Gstar]=ML_edynfactors2(ustar,q);                                    % Estimating Dynamic Factors Innovation
    for ii=1:h; CLstar(:,:,ii)=(lambda*Bstar(:,:,ii)*Gstar).*WWW; end;      % MA Representation of common components
    BB = CumImp(CLstar, cd);                                                % Comulates IRF if necessary    
    [H eflag(bb,1)]=ML_SLRrestrictions(ID,BB,BB(:,:,end),K0,K1);            % Rotation Matrix    
    CC(:,:,:,bb) = ML_ComputeIrf(BB,H);                                     % Compute IRF           
end
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

