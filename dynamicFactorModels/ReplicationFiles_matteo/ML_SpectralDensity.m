%% ML_SpectralDensity - Computes the spectral density matrix
%
% Sigma = ML_SpectralDensity(X, m, h)
%

% Written by Matteo Luciani (matteo.luciani@ulb.ac.be)
% Modified version of spectral.m by Matteo Barigozzi

function [Sigma, omega] = ML_SpectralDensity(X, m, h)

%% Preliminary settings
[T,n] = size(X);

if nargin==1; m = floor(sqrt(T)); h = m; end
if nargin==2; h = m; end

%% Compute M covariances
M = 2*m+1;
B = triang(M);                                                              % Triangular window (similar Bartlett)
Gamma = zeros(n,n,M);
for k = 1:m+1,
    Gamma(:,:,m+k) = B(m+k)*(X(k:T,:))'*(X(1:T+1-k,:))/(T-k);
    Gamma(:,:,m-k+2) = Gamma(:,:,m+k)';
end

%% Compute the spectral density matrix in H points
H = 2*h+1;
omega=-2*pi*h/H:2*pi/H:2*pi*h/H;
Factor = exp(-sqrt(-1)*(-m:m)'*(omega));
Sigma = zeros(n,n,H);
for j = 1:n
    Sigma(j,:,:) = squeeze(Gamma(j,:,:))*Factor;
end