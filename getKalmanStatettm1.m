function [attm1,pttm1]=getKalmanStatettm1(ptt,att,H,covErr)
% Do not get confused with notation (attm1 and pttm1 are the prediction for the next
% time step based on current estimate att, ptt, and the state transition matrix H)
%
% The state space representation format
%   Z_t=eigVec*a_t+\nu+\eps
%   a_t=(eigVec'*B)*a_(t-1)+eigVec'*\eta
%   a_t=H*a_(t-1)+eigVec'*\eta
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


attm1=H*att;
pttm1=H*ptt*H'+covErr;

end