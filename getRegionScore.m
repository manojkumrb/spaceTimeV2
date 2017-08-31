function [selIndex]=getRegionScore(tol, devRegion, varRegion)
% Calculates the probability of a region being out of tolerance and finds
% the entropy (average information) obtained by measuring a given region
%
%
% The predicted deviations are assumed to be from a gaussian distribution
% with a mean at the predicted value and variance equal to that found by
% the kalman filter
%
% tol       = the defined tolerance vlaue
% devRegion = the predicted or filterd value of the deviation
% varRegion = The variance of the predicted or filterd value of the deviation
% selIndex  = entropy gain obtained be measuring the region

if length(devRegion)~=length(varRegion)
    display('imput deviations and corresponding variances length donot match');   
    return;    
end

probOut=1-cdf('Normal',tol,devRegion,sqrt(varRegion))+cdf('Normal',-tol,devRegion,sqrt(varRegion));
probOut(probOut<=eps)=1; % to avoid the infinity case as prob becomes zero and make the it equal to the limiting value as p-->0
selIndex=sum(-probOut.*log(probOut))/length(devRegion);  % normalised to counter the effect of number of points in the region


end