function [ytt,VarYtt,ytp1,VarYtp1]=getPredictionsNew(eigVec,att,ptt,H,V,R,covErr,hyp,nodeCoordinates,xTemp,zTemp,tempEigVec)
% Includes krigging the error as i had skipped it in the first code
% gives filtered and predicted values for deviation and its variance 
%
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
% eigVec    = is the system matrix \phi .. which is equivalent to eigen
%             vector in our case
% ptt       = state variable covariance matrix (of size k x k) at t given (t)
% att       = state variable vector of length k, at at t given (t)
% ytt       = the filtered value for deviation, (in the code krigging part has not been taken into account)
% VarYtt    = variance for the filterd value
% yttp1     = the Predicted value for deviation in (t+1)
% VarYtt    = variance for the predicted value in (t+1)
% *Temp		= all data related to measured nodes


% simple krigging the residuals
m0 = {'meanZero'};  hyp.mean = [];      % no hyperparameters are needed
covFunc = {@covMaternard,3};  % Matern class d=3
likfunc = {'likGauss'};  % gaussian likelihood


yTemp=att'*tempEigVec';
diff=zTemp-yTemp;
[mpredicted,varPred] = gp(hyp, @infExact, m0, covFunc, likfunc, xTemp,diff', nodeCoordinates );

ytt=att'*eigVec'+mpredicted';   % taking transpose to maintain convention of deviation matrix(i.e. rep x nodes)
VarYtt=diag(eigVec*ptt*eigVec'+(V+R))+varPred;
ytp1=att'*H'*eigVec';  % taking transpose to maintain convention of deviation matrix
VarYtp1=diag(eigVec*(H*ptt*H'+covErr)*eigVec'+(V+R));

end