% BLL_NumberFactorsShocks - Determines the Number of Factors and Shocks
%

clear all; close all; clc; 
ML_graph_options
[Data, Label, Name, Dates, cd] = BLL_ReadData(2,1,1);                       % Load data
x=ML_Standardize(Data);                                                     % standardize data
[T, N]=size(x);                                                             % size of the dataset
kmax=10; cmax=3; npace=1; m = floor(.5*T^(1/2)); nbck=floor(N/3);           % parameters for IC on number of shocks

ML_nfactorsBN(x,20);                                                        % Bai and Ng (2002) IC for number of Factors
Sigma=ML_SpectralDensity(x,m,m); Sigma(:,:,1:7)=[];
for ii=1:8; autoval(:,ii)=real(eig(Sigma(:,:,ii))); end

disp('   *****************')
disp(['    Eigenvalues ' char(931) '(' char(952)  ')'])
disp('   -----------------')
disp(['      ' char(955) '(0)      '    char(955) '(' char(952)  ')' ])
disp([autoval(1:20,1)  mean(autoval(1:20,:),2)])
disp('   *****************')

nrep=100;                                                                   % number of times we want to run the criteria for n shocks
[qhat1, qhat2]=BLL_Criteria2(x, m, kmax, nbck, cmax, nrep);
