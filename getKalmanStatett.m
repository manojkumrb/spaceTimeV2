function [kt,att,ptt]=getKalmanStatett(pttm1,attm1,eigVec,R,V,zt)
% 
% The state space representation format
%   Z_t=eigVec*a_t+\nu+\eps
%   a_t=(eigVec'*B)*a_(t-1)+eigVec'*\eta
%   a_t=H*a_(t-1)+eigVec'*\eta
%
% zt        = measurement vector of length n, at time t
% pttm1     = state variable covariance matrix (of size k x k) at t given (t-1)
% attm1     = state variable vector of length k, at at t given (t-1)
% R         = covariance matrix of the \eps process
% V         = covariance matrix of the \nu process
% kt        = kalman gain at time t
% ptt       = state variable covariance matrix (of size k x k) at t given (t)
% att       = state variable vector of length k, at at t given (t)


kt=(pttm1*eigVec')/(R+V+eigVec*pttm1*eigVec');
att=attm1+kt*(zt-eigVec*attm1);
ptt=pttm1-kt*eigVec*pttm1;


end