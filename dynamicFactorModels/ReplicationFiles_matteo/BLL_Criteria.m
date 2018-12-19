% BLL_Criteria - Determines the Number of Common Shocks and Common Trends
% 
%  This function implement the information criteria described in Barigozzi,
%  Lippi and Luciani (2016) to determin the number of common trends in
%  large panels. 
%  This function also implement the information criteria of 
%  Hallin and Liska (2007 - JASA)
% 
% [qhat1, qhat2]=BLL_Criteria(X, M, kmax, nbck, cmax, penalty,)
% 
% Inputs:
% 	   X - (T x n) stationary data matrix
%   kmax - maximum number of shocks
%      M - covariogram truncation, e.g. M = 0.75*sqrt(T); M = 0.5*sqrt(T);
%   nbck - number of sublocks to be used
% Outputs:
%  qhat1 - determines the number of shocks using a large window
%  qhat2 - determines the number of shocks using a small window
% 
% Written by Matteo Luciani (matteo.luciani@frb.gov)
% 
% Reference: Barigozzi, Matteo, Marco Lippi, and Matteo Luciani (2016).   
%       “Non-Stationary Dynamic Factor Models for Large Datasets,”  
%       Finance and Economics Discussion Series 2016-024.   
% 

function [qhat1, qhat2]=BLL_Criteria(X, M, kmax, nbck, cmax)

npace=1; 
step=500; 

[T,n] = size(X);                                                            % Size of the datatset
x = ML_Standardize(X);                                                      % Standardize
Eigenvalues=NaN(length(n-nbck:npace:n),kmax,2);                             % Preallocates
a1=kmax+1; a2=2;

    
% Running the Criteria %%
s=0;
for N = n-nbck:npace:n
    s=s+1; 
    [~, Ns]=sort(rand(n,1));
    eigv = EigenvaluesSpectralDensity(x(1:T,Ns(1:N)),M);                    % Eigenvalues of the Spectra Density Matrix        
    IC1 = flipud(cumsum(flipud(eigv))); IC1 = IC1(1:kmax+1,:);              % cumulative sum of eigenvalues in decreasing order    
    for jj=1:a2;Eigenvalues(N-n+nbck+1,:,jj)=eigv(1:10,jj)';end             % Store the eigenvalues
    
    for pp=1:4;                                                             % computes the criterion using different penalties
        p=SelectPenalty(pp,T,M,N,a1,a2);
        T0=repmat((0:kmax)',1,a2).*p;
        for c = 1:floor(cmax*step);            
            cc = c/step;
            IC_log = log(IC1./N) + T0*cc;   [~, rr(1,:)]=min(IC_log);       % Log criterion                                    
            for jj=1:2; bll{1,pp}(s,c,jj)=rr(1,jj)-1; end                  % -------------
        end
    end    
end 

%%% ----------------------------------------- %%%
%%% Select Automatically the Number of Shocks %%%
%%% ----------------------------------------- %%%

cr=(1:floor(cmax*500))'/500;

for pp=1:4;    
    for jj=1:2;
        BLL{1,jj,pp}(1,1)=kmax;
        BLL{1,jj,pp}(1,2:3)=0;
    end;
    kk=1;
    
    sbll{1,pp}=squeeze(std(bll{1,pp}));
    for jj=1:2;                                                             % common shocks or common trends
        c1=2;
        for ii=1:size(cr,1);
            if sbll{1,pp}(ii,jj)==0;                                          % If the number of shocks is always the same across blocks
                if bll{1,pp}(end,ii,jj)==BLL{1,jj,pp}(c1-1,1);
                    BLL{1,jj,pp}(c1-1,3)=cr(ii);
                else
                    BLL{1,jj,pp}(c1,1)=bll{1,pp}(end,ii,jj);
                    BLL{1,jj,pp}(c1,2:3)=cr(ii);
                    c1=c1+1;
                end;
            end;
        end;
        BLL{1,jj,pp}(:,4)=BLL{1,jj,pp}(:,3)-BLL{1,jj,pp}(:,2);              % Computes Window Size
        q=BLL{1,jj,pp}(find(BLL{1,jj,pp}(2:end,4)>.1)+1,1);                 % Number of Shocks with Large Window
        qhat1(pp,kk)=q(1);                                                  % ----------------------------------
        q=BLL{1,jj,pp}(find(BLL{1,jj,pp}(2:end,4)>.01)+1,1);                % Number of Shocks with Small Window
        qhat2(pp,kk)=q(1);                                                  % ----------------------------------
        kk=kk+1;
    end    
end


% Display Results %%

T1={'BLL log','HL log'};
T2=['a' 'b' 'c' 'd'];
for jj=[1 3:5 7:9]; ZZ{jj,1}=' '; end; ZZ{1,2}='Penalty ';
for jj=1:4; ZZ{jj+1,2}=T2(jj); ZZ{jj+5,2}=T2(jj); end
ZZ{2,1}='Large Window';
ZZ{6,1}='Small Window';
Table=[];
for jj=1:2; Table=cat(2,Table,repmat('   ',9,1),char(ZZ{:,jj})); end

for kk=1:2;
    ZZ{1,kk+2}=T1{kk};
    for pp=1:4;
        ZZ{1+pp,kk+2}=num2str(qhat1(pp,kk));
        ZZ{5+pp,kk+2}=num2str(qhat2(pp,kk));
    end
    Table=cat(2,Table,repmat('   ',9,1),char(ZZ{:,kk+2}));
end

T3=size(Table,2);
disp('')
disp(repmat('*',1,T3))
disp('       Number of Shocks Selected with Automatic Procedure')
disp(repmat('-',1,T3))
disp(Table(1,:));
disp(repmat('-',1,T3))
disp(Table(2:5,:));
disp(repmat('-',1,T3))
disp(Table(6:9,:));
disp(repmat('*',1,T3))

disp('a = ((M/T)^0.5 + M^(-2) + N^(-1))*log(min([(T/M)^0.5;  M^2; N]))')
disp('b = (min([(T/M)^0.5;  M^2; N])).^(-1/2)')
disp('c = (min([(T/M)^0.5;  M^2; N])).^(-1)*log(min([(T/M)^0.5;  M^2; N]))')
disp('d = (min([(T/M)^0.25;  M^2; N])).^(-1)*log(min([(T/M)^0.25;  M^2; N]))')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EigenvaluesSpectralDensity - Computes Eigenvalues of the spectral density matrix

function E = EigenvaluesSpectralDensity(x,M)

w=M; 
Sigma = ML_SpectralDensity(x,M); Sigma(:,:,1:M)=[];                         % Estimate spectral density matrix at all frequencies
for jj=1:size(Sigma,3); eigenvalues(:,jj)=real(eig(Sigma(:,:,jj))); end;    % Eigenvalues at different frequencise
E(:,1) = eigenvalues(:,1);                                                  % Eigenvalues at the zero frequency
E(:,2) = [eigenvalues(:,1)  eigenvalues(:,2:w+1)*2]*ones(w+1,1)/(2*w+1);    % eigenvalues of spectral density matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ML_SpectralDensity - Computes the spectral density matrix

function Sigma = ML_SpectralDensity(X, m, h)

% Preliminary settings
[T,n] = size(X);

if nargin==1; m = floor(sqrt(T)); h = m; end
if nargin==2; h = m; end

% Compute M covariances
M = 2*m+1;
B = triang(M);                                                              % Triangular window (similar Bartlett)
Gamma = zeros(n,n,M);
for k = 1:m+1,
    Gamma(:,:,m+k) = B(m+k)*(X(k:T,:))'*(X(1:T+1-k,:))/(T-k);
    Gamma(:,:,m-k+2) = Gamma(:,:,m+k)';
end

% Compute the spectral density matrix in H points
H = 2*h+1;
Factor = exp(-sqrt(-1)*(-m:m)'*(-2*pi*h/H:2*pi/H:2*pi*h/H));


Sigma = zeros(n,n,H);
for j = 1:n; Sigma(j,:,:) = squeeze(Gamma(j,:,:))*Factor; end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y, M, s] = ML_Standardize(x)

T=size(x,1);
s = std(x);
M = mean(x);
ss = ones(T,1)*s;
MM = ones(T,1)*M;
y = (x-MM)./ss;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p=SelectPenalty(J,T,M,N,a1,a2)

if J == 1;
    p = ((M/T)^0.5 + M^(-2) + N^(-1))*log(min([(T/M)^0.5;  M^2; N]))*ones(a1,a2);
elseif J == 2;
    p = (min([(T/M)^0.5;  M^2; N])).^(-1/2)*ones(a1,a2);
elseif J == 3;
    p = (min([(T/M)^0.5;  M^2; N])).^(-1)*log(min([(T/M)^0.5;  M^2; N]))*ones(a1,a2);
elseif J == 4;    
    p = (min([(T/M)^0.25;  M^2; N])).^(-1)*log(min([(T/M)^0.25;  M^2; N]))*ones(a1,a2);
end

