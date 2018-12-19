% BLL_Criteria2 - Determines the Number of Common Shocks and Common Trends
% 
%  This function implements BLL_Criteria nrep times
% 
% [qhat1, qhat2]=BLL_Criteria(X, M, kmax, nbck, cmax, penalty,nrep)
% 
% Inputs:
% 	   X - (T x n) stationary data matrix
%   kmax - maximum number of shocks
%      M - covariogram truncation, e.g. M = 0.75*sqrt(T); M = 0.5*sqrt(T);
%   nbck - number of sublocks to be used
%   nrep - number of replications
% 
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

function [qhat1, qhat2]=BLL_Criteria2(X, M, kmax, nbck, cmax,nrep)

tic

% Running the Criteria %%
matlabpool size;
if ans==0; matlabpool; end

t1=NaN(4,2,nrep); t2=t1;

parfor nn=1:nrep;        
    [t1(:,:,nn), t2(:,:,nn)]=BLL_Criteria0(X, M, kmax, nbck, cmax);
end

for nn=1:nrep;
    for pp=1:4;
        for kk=1:2;    % BLL & HL         
            qhat1(nn,pp,kk)=t1(pp,kk,nn);
            qhat2(nn,pp,kk)=t2(pp,kk,nn);
        end
    end
end

T1={'BLL log','HL log'};
T2={'q','p1','p2','p3','p4'};

for kk=1:2; 
    Table=[];
    Z1=ML_tabulate2(qhat1(:,:,kk));
    for jj=1:5;
        ZZ{1,jj}=T2{jj};
        for ii=1:size(Z1,1);
            if jj==1; ZZ{ii+1,jj}=num2str(Z1(ii,jj),'%7.0f'); 
            else ZZ{ii+1,jj}=num2str(Z1(ii,jj),'%7.2f'); end
        end
        Table=cat(2,Table,repmat('  ',size(ZZ,1),1),char(ZZ(:,jj)));
    end
    disp(' ')
    disp(repmat('*',1,size(Table,2)))
    disp(['  Large Window - ' T1{kk}])
    disp(Table(1,:))    
    disp(repmat('-',1,size(Table,2)))
    disp(Table(2:end,:));
    disp(repmat('-',1,size(Table,2)))
    clear ZZ
end

for kk=1:2;
  Table=[];
    Z1=ML_tabulate2(qhat2(:,:,kk));
    for jj=1:5;
        ZZ{1,jj}=T2{jj};
        for ii=1:size(Z1,1);
            if jj==1; ZZ{ii+1,jj}=num2str(Z1(ii,jj),'%7.0f'); 
            else ZZ{ii+1,jj}=num2str(Z1(ii,jj),'%7.2f'); end
        end
        Table=cat(2,Table,repmat('  ',size(ZZ,1),1),char(ZZ(:,jj)));
    end
    disp(' ')
    disp(repmat('*',1,size(Table,2)))
    disp(['   Small Window - ' T1{kk}])
    disp(Table(1,:))    
    disp(repmat('-',1,size(Table,2)))
    disp(Table(2:end,:));
    disp(repmat('-',1,size(Table,2)))
    clear ZZ
end



disp('p1 = ((M/T)^0.5 + M^(-2) + N^(-1))*log(min([(T/M)^0.5;  M^2; N]))')
disp('p2 = (min([(T/M)^0.5;  M^2; N])).^(-1/2)')
disp('p3 = (min([(T/M)^0.5;  M^2; N])).^(-1)*log(min([(T/M)^0.5;  M^2; N]))')
disp('p4 = (min([(T/M)^0.25;  M^2; N])).^(-1)*log(min([(T/M)^0.25;  M^2; N]))')


toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [qhat1, qhat2]=BLL_Criteria0(X, M, kmax, nbck, cmax)
s=0;    
npace=1;
step=500;

[T,n] = size(X);                                                            % Size of the datatset
x = ML_Standardize(X);                                                      % Standardize
a1=kmax+1; a2=2;

for N = n-nbck:npace:n
    s=s+1;
    [~, Ns]=sort(rand(n,1));
    eigv = EigenvaluesSpectralDensity(x(1:T,Ns(1:N)),M);                    % Eigenvalues of the Spectra Density Matrix
    IC1 = flipud(cumsum(flipud(eigv))); IC1 = IC1(1:kmax+1,:);              % cumulative sum of eigenvalues in decreasing order
    for pp=1:4;                                                             % computes the criterion using different penalties
        p=SelectPenalty(pp,T,M,N,a1,a2);
        T0=repmat((0:kmax)',1,a2).*p;
        for c = 1:floor(cmax*step);
            cc = c/step;
            IC_log = log(IC1./N) + T0*cc;   [~, rr(1,:)]=min(IC_log);       % Log criterion
            for ll=1;
                for jj=1:2;
                    bll{ll,pp}(s,c,jj)=rr(ll,jj)-1;                         % -------------
                end
            end
        end
    end
end

%%% ----------------------------------------- %%%
%%% Select Automatically the Number of Shocks %%%
%%% ----------------------------------------- %%%

cr=(1:floor(cmax*500))'/500;

for pp=1:4;
    for ll=1;
        for jj=1:2;
            BLL{ll,jj,pp}(1,1)=kmax;
            BLL{ll,jj,pp}(1,2:3)=0;
        end;
    end
    
    kk=1;
    for ll=1;                                                               % common shocks or common trends
        sbll{ll,pp}=squeeze(std(bll{ll,pp}));
        for jj=1:2;                                                         % log or no-log
            c1=2;
            for ii=1:size(cr,1);
                if sbll{ll,pp}(ii,jj)==0;                                   % If the number of shocks is always the same across blocks
                    if bll{ll,pp}(end,ii,jj)==BLL{ll,jj,pp}(c1-1,1);
                        BLL{ll,jj,pp}(c1-1,3)=cr(ii);
                    else
                        BLL{ll,jj,pp}(c1,1)=bll{ll,pp}(end,ii,jj);
                        BLL{ll,jj,pp}(c1,2:3)=cr(ii);
                        c1=c1+1;
                    end;
                end;
            end;
            BLL{ll,jj,pp}(:,4)=BLL{ll,jj,pp}(:,3)-BLL{ll,jj,pp}(:,2);                    % Computes Window Size
            q=BLL{ll,jj,pp}(find(BLL{ll,jj,pp}(2:end,4)>.1)+1,1);                     % Number of Shocks with Large Window
            qhat1(pp,kk)=q(1);                                                   % ----------------------------------
            q=BLL{ll,jj,pp}(find(BLL{ll,jj,pp}(2:end,4)>.01)+1,1);                    % Number of Shocks with Small Window
            qhat2(pp,kk)=q(1);                                                   % ----------------------------------
            kk=kk+1;
        end
    end
end



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

% 
% disp('')
% disp('***********************************************')
% disp('**    Barigozzi, Lippi & Luciani Criteria    **')
% disp('**       for number of common trends         **')
% disp('***********************************************')
% disp('')
% disp('             No-log criterium')
% disp('  n trends   start        end      size')
% disp(BLL{1,1,2})
% disp('')
% disp('             log criterium')
% disp('    trends   start        end      size')
% disp(BLL{1,2,2})
% 
% disp('')
% disp('************************************************')
% disp('**                                            **')
% disp('**    Hallin & Liska (2007, JASA) Criteria    **')
% disp('**                                            **')
% disp('************************************************')
% disp('')
% disp('             No-log criterium')
% disp('    shocks   start        end      size')
% disp(BLL{2,1,2})
% disp('')
% disp('             log criterium')
% disp('    shocks   start        end      size')
% disp(BLL{2,2,2})
