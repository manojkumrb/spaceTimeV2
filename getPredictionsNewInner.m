function [ytt,VarYtt,ytp1,VarYtp1,yttU,ytp1U,VarYttU]=getPredictionsNewInner(interpEigVec,...
						eigVeckey,tempEigVec,att,ptt,H,V,R,covErr,hyp,nodeCoordinatesU,...
						nodeCoordinatesM,xTemp,zTemp)
% Includes krigging the error as i had skipped it in the first code
% gives filtered and predicted values for deviation and its variance 
%
%
% The state space representation format
%   Z_t=eigVec*a_t+\nu+\eps
%   a_t=(eigVec'*B)*a_(t-1)+eigVec'*\eta
%   a_t=H*a_(t-1)+eigVec'*\eta
%
% zt				= measurement vector of length n, at time t
% pttm1				= state variable covariance matrix (of size k x k) at t given (t-1)
% attm1				= state variable vector of length k, at at t given (t-1)
%
% interpEigVec		= eigen vectors corresponding to all the points in the mesh
% eigVecKey			= is the system matrix \phi .. which is equivalent to eigen
%					  vectors corresponding to the **key points** in our case
% *Temp,temp*		= all data related to measured nodes (partial measurements)
% att				= state variable vector of length k, at at t given (t)
% ptt				= state variable covariance matrix (of size k x k) at t given (t)
% H					= the state transition matrix
% R					= covariance matrix of the \eps process
% V					= covariance matrix of the \nu process
% Q					= covariance matrix of the \eta process
% covError			= the state transition equation error covariance matrix (eigVec'*Q*eigVec)
% hyp				= optimised hyperparamters for krigging
% nodeCoordinatesU	= coordinates of all mesh points
% nodeCoordinatesM	= coordinates of all measured mesh points, corresponds
%					  to key points
%
%
% ytt				= the filtered value for deviation for the key points
% VarYtt			= variance for the filterd value
% yttp1				= the Predicted value for deviation in (t+1) for key points
% VarYtt			= variance for the predicted value in (t+1)
% yttU				= the filtered value for deviation of all mesh points
% yttp1				= the Predicted value for deviation in (t+1)of all mesh points



% simple krigging the residuals
m0 = {'meanZero'};  hyp.mean = [];      % no hyperparameters are needed
covFunc = {@covMaternard,5};  % Matern class d=3
likfunc = {'likGauss'};  % gaussian likelihood


yTemp=att'*tempEigVec';
diff=zTemp-yTemp;
[mpredicted,varPred] = gp(hyp, @infExact, m0, covFunc, likfunc, xTemp,diff', nodeCoordinatesM );

ytt=att'*eigVeckey'+mpredicted';   % taking transpose to maintain convention of deviation matrix(i.e. rep x nodes)
VarYtt=diag(eigVeckey*ptt*eigVeckey'+(V+R))+varPred;
ytp1=att'*H'*eigVeckey';  % taking transpose to maintain convention of deviation matrix
VarYtp1=diag(eigVeckey*(H*ptt*H'+covErr)*eigVeckey'+(V+R));

% prediction for all points on the mesh 
% variance is neglected for now.. will consider it later.. variance for
% only key points is considered now.
[mpredictedU,varPredU] = gp(hyp, @infExact, m0, covFunc, likfunc, xTemp,diff', nodeCoordinatesU );
yttU=att'*interpEigVec'+mpredictedU';   % taking transpose to maintain convention of deviation matrix(i.e. rep x nodes)
VarYttU=sum(interpEigVec.*(ptt*interpEigVec')',2)+varPredU; % computation to calculate diagonal values of product of 3 matrices
% this variance neglects other kalman filter variance sources (V + R) as in line 54
ytp1U=att'*H'*interpEigVec';  % taking transpose to maintain convention of deviation matrix
% VarYtp1U=diag(interpEigVec*(H*ptt*H')*interpEigVec');


end