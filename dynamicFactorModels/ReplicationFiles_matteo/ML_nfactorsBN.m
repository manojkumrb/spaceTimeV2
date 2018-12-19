% ML_nfactorsBN - Runs, and shows results of  Bai & Ng (2002) criteria

% Written by Matteo Luciani (matteo.luciani@ulb.ac.be)

function [PR IC PC]=ML_nfactorsBN(y,rmax)

for jj=1:4; [a b PR(:,jj) IC(:,jj) PC(:,jj)]=ML_nfactors(y,rmax,jj); close; R(jj,:) = [a b]; end;   % Bai&Ng (2002) Criteria for Static Factors
disp('-------------------------------------------------------------');
disp('Selected Number of Static Factors with Bai & Ng 2002 Criteria');
disp('    IC    PC');
disp(R);
disp('-------------------------------------------------------------');
clear a b c R;

% ML_nfactors - Determine the Number of Static Factors as Bai and Ng (2002, Econometrica)
%
% Computes IC, PC & BIC Criteria
%
% [ICx, PCx]=ML_NFACTORS_b(y,kmax,jj);
% Inputs:
%   y is observed (data matrix)
%   kmax = maximum nuber of factors to test
%   jj = criteria {1,2,3,4}
% Outputs:
%   ICx  = Selected Number of Static Factors with the IC criteria
%   PCx  = Selected Number of Static Factors with the PC criteria
%

% Written by Matteo Luciani (matteo.luciani@ulb.ac.be)
% This is a customized version of the codes available on Serena Ng webpage

function [ICx, PCx, PR, IC, PC]=ML_nfactors(y,kmax,jj)

[T,N]=size(y);
NT=N*T; NT1=N+T; GCT=min([N;T]);
ii=1:1:kmax; id=1:N;
CT1=zeros(1,kmax);                                                          % Penality Function for IC
CT2=zeros(1,kmax);                                                          % Penality Function for PC
Sigma=zeros(1,kmax+1);                                                      % V(k,F)
IC=zeros(size(CT1,1),kmax+1);                                               % Panel Information Criteria
PC=zeros(size(CT2,1),kmax+1);                                               % Panel Cp
RV=zeros(size(CT1,1),kmax);                                                 % Percentage of Total Variance
PR=zeros(size(CT1,1),kmax);                                                 % Proportions of Variance

% ---------------------------------------------- %
% Factors Estimation and Computation of Criteria %
% ---------------------------------------------- %
w=((T-1)^(-1));                                                             % This is just a shorcut for computing variances
[Fhat0,eigval]=svd(w*y'*y);

%%% ======================= %%%
%%% Computation of Criteria %%%
%%% ======================= %%%

lambda=Fhat0(:,1:kmax)*sqrt(N); fhat=y*lambda/N; ehat=y-fhat*lambda';
Sigma(kmax)=mean(sum(ehat.*ehat/T));                                        % sigma hat square
for i=1:1:kmax;
	lambda=Fhat0(:,1:i)*sqrt(N); fhat=y*lambda/N; ehat=y-fhat*lambda';
	Sigma(i)=mean(sum(ehat.*ehat/T));
	y_sigma=w*y'*y;                                                                                    % Variance-Covariance Matrix of the Differenced Data
	if     jj==1;  CT1(1,i)=i*(NT1/NT)*log(NT/NT1); CT2(1,i)=i*Sigma(kmax)*(NT1/NT)*log(NT/NT1);       % Computation of Penalities for IC
	elseif jj==2;  CT1(1,i)=i*(NT1/NT)*log(GCT);    CT2(1,i)=i*Sigma(kmax)*(NT1/NT)*log(GCT);          % --------------------------------
	elseif jj==3;  CT1(1,i)=i*(log(GCT)/(GCT));     CT2(1,i)=i*Sigma(kmax)*(log(GCT)/(GCT));           % --------------------------------
	else           CT1(1,i)=i*(N+T-i)*log(NT)/NT;   CT2(1,i)=i*Sigma(kmax)*(N+T-i)*log(NT)/NT;    end; % --------------------------------
	IC(:,i)=log(Sigma(i))+CT1(:,i);
	PC(:,i)=Sigma(i)+CT2(:,i);
	PR(:,i)=eigval(i,i)/trace(y_sigma);
	RV(:,i)=trace(eigval(1:i,1:i))/trace(y_sigma);
end;

% ----------------- %
% Output Displaying %
% ----------------- %

eigval_g=diag(eigval); ee=1:kmax;
figure,plot(ee,eigval_g(1:kmax),'k-*'), axis tight
title('Value of the Eigenvalues','FontWeight', 'bold', 'FontName', 'times'),
xlabel('N° of Eigenvalues', 'FontName', 'times');
set(0,'DefaultLineLineWidth',1);
Sigma(kmax+1)=mean(sum(y.*y/T));
IC(:,kmax+1)=log(Sigma(kmax+1));
PC(:,kmax+1)=Sigma(kmax+1);

if jj==1;
	disp('**************************************************');
	disp('**                                              **');
	disp('**   Determining the Number of Static Factors   **');
	disp('**                                              **');
	disp('**************************************************');
	fprintf('N = ');
	disp(N);
	fprintf('T = ');
	disp(T);
	fprintf('rmax = ');
	disp(kmax);
	%fprintf('N° Factors ');
	%disp([(1:1:kmax) 0]);
end;

outputdisplay=0;
if outputdisplay==1;
	if      jj==1; disp('       IC1       PC1');
	elseif  jj==2; disp('       IC2       PC2');
	elseif  jj==3; disp('       IC3       PC3');
	else           disp('       IC4       PC4');
	end;
	disp([IC' PC']);
end;

if jj==4;
	for ii=1:kmax;
		ar(ii,1)=ML_autocorrelation(fhat(:,ii),1);
	end;
	disp('N° Factors    M.R-sq    A.R-sq    AR(1)');
	disp([(1:kmax)' PR' RV' ar]);
end;

%%% ------------------------------------------- %%%
%%% Selected Number of Factors by each Criteria %%%
%%% ------------------------------------------- %%%
ICx=0; PCx=0;
for i=1:kmax;
	if IC(:,i)==min(IC); ICx=i;  end;
	if PC(:,i)==min(PC); PCx=i;  end;
end;