% BLL_SLRImp - Computes IRF for I(1) DFM with SLR restrictions (VECM)

% Written by Matteo Luciani (matteo.luciani@ulb.ac.be)

function [C,e] = BLL_SLRImp(Y,cd,q,d,r,p,h,ID,K0,K1)

y=ML_diff(Y);                                                               % 1st Differences
[yy, ~, sy]=ML_Standardize(y);                                              % Standardize
X=ML_detrend(Y,1);                                                          % Detrend Data
Z=X./repmat(sy,size(X,1),1);                                                % BLL - Standardization 
[T, N]=size(y);                                                             % Panel Dimensions
[~,lambda]=ML_efactors2(yy,r,2);                                            % estimate loadings
F=Z*lambda/N;                                                               % BLL Factors
[~, ~, ~, CL, ~, v]=BLL_VECM(F,p,r-q+d,h,1);                                % Estimate the VECM
[V D] = eigs(cov(v),q,'lm'); G=V*(D^.5); u=v*V*(D^-.5);                     % common shocks
for jj=1:h; Phi(:,:,jj)=(lambda*CL(:,:,jj)*G).*repmat(sy',1,q); end;        % Impulse Responses
Phi = CumImp(Phi,cd);                                                       % comulate IRFs if it is necessary
H=ML_SLRrestrictions(ID,Phi,Phi(:,:,end),K0,K1);                            % Rotation Matrix
for jj=1:h; C(:,:,jj)=Phi(:,:,jj)*H; end;                                   % Structural IRF
e=u*H;                                                                      % StructuralShocks


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






