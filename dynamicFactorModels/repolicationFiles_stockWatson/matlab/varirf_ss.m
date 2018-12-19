function irf = varirf_ss(coef,i,hor)

% MODEL: Q, M, G for model written in SS form
%           y(t) = Q*z(t)
%           z(t) = M*z(t-1) + G*u(t)
%           var(u(t)) = I
% coef: coefficients {Q, M, G}.
% i   : vector with shocks of interest
% hor : horizon for impulse response
% 
% -- OUTPUT --
% irf(i,j,t): impulse response of variable i to shock j at horizon t-1.

% PRELIMINARIES
Q = coef.Q;    M = coef.M;    G = coef.G;    % coefficients
ndim = size(Q,1);                               % data dimensions
if i == 0        
    nshocks = size(G,2);	% number of shocks to compute impulse response for
    i = 1:nshocks;          % compute irf for all shocks
else
    nshocks = length(i);
end

% IMPULSE RESPONSES
irf = NaN(ndim,nshocks,hor);
for j = 1:nshocks;
    i_shock = i(j);
    x = G(:,i_shock);
    for i_hor = 1:hor;
        y = Q*x;
        irf(:,j,i_hor) = y;
        x = M*x;
    end;
end;

end