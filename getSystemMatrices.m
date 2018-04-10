function [eigVec,R,V,H,covErr]=getSystemMatrices(nodeCoord,dev,sigmaMes,nu)


%% the kalman filter notations are from Cressie and wikle (biometrica)
%
% nodeCoord = nx 3 matrix of point coordinates in teh same order as defined
%             in the domain of the part
% dev       = t x n matrix of series of deviations to be modelled. The deviations
%             are centered. t id the number of time instances, n is the
%             number of points
% sigmaMes  = measuremetn error standard deviation
% nu        = a structure representing small scale spatial variation which
%             donot have a temprally dynamic structure
% nu.Type   = 'diag' or 'covFunc'
% nu.CovPar = parameters of the cov function if the type is covFunc, the
%             function by default is maternARD with d=3, needs GP toolbox
% nu.StdDev = standard deviation value if the type is diag
%
% State equations
%       Z_t=eigVec*a_t+\nu+\eps
%       a_t=(eigVec'*B)*a_(t-1)+eigVec'*\eta
%          =H*a_(t-1)+eigVec'*\eta
%
% zt        = measurement vector of length n, at time t
% pttm1     = state variable covariance matrix (of size k x k) at t given (t-1)
% attm1     = state variable vector of length k, at at t given (t-1)
% H         = the state transition matrix
% R         = covariance matrix of the \eps process
% V         = covariance matrix of the \nu process
% Q         = covariance matrix of the \eta process
% covError  = the state transition equation error covariance matrix (eigVec'*Q*eigVec)
% ptt       = state variable covariance matrix (of size k x k) at t given (t)
% att       = state variable vector of length k, at at t given (t)
%
%dependencies: GP toolbox, nearestspd,



% to obtain major deviation patterns
[~ ,singVal,eigVec]=svds(dev,25);
% eigVal=diag(singVal).^(2)./size(dev,2);
% pc=dev*eigVec;%plot(pc);                    % principal components
% x=pc*eigVec';
% diff=dev-x;

% estimation of covariances c_z
[n, t]=size(dev);
% zero lag covariance
gMean=sum(sum(dev))/(n*t);
C0z     =((dev-gMean)'*(dev-gMean))./t;
% lag one covariance
devCenterdT=dev(2:end,:);
devCenterdTminus1=dev(1:end-1,:);
C1z=((devCenterdT-gMean)'*(devCenterdTminus1-gMean))./(t-1);
C1z = nearestSPD(C1z);                          % finding nearest positive definite matrix

% measurement error variance
R=sigmaMes.*eye(size(C1z));

% coloured error process covariance matrix
if strcmpi(nu.Type,'diag')
    V=eye(size(R))*nu.StdDev;
elseif strcmpi(nu.Type,'covFunc')
    covFunc = {@covMaternard,3};   % Matern class d=3
    V=feval(covFunc{:},nu.CovPar,nodeCoord); % optimised values obtained from average of maximum likelihood values of all samples
    % V=feval(@covMaternard,3,hyp.cov,nodeCoord); % this method to evaulate cov func also works
else
    display('unknown noise process defined')
    return;
end


% zero lag y process covariance
C0y=C0z-R;
C0y = nearestSPD(C0y);

% B matrix ie. the matrix that specifies autocorrelation b/w a_t and a_(t-1)
J=eigVec'; % I take eigVec'eigVec to be identity matrix as eigVec is orthogonal (eigVec corresponds to phi the eigen vector matrix in teh source paper)
B=(C1z*J')/(J*(C0z-R-V)*J');
% estimation of JQJ' matrix ie. the state error covariance matrix
H=J*B;
covErr=J*(C0z-R-V)*J'-J*C1z*J'*H'-H*J*C1z'*J'+H*J*(C0z-R-V)*J'*H';
covErr = nearestSPD(covErr);



end