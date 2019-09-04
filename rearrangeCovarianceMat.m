function covOut=rearrangeCovarianceMat(covOrig,nodeIDoutRect)
% covOrig       = n x n covariance matrix
% nodeIDOutRect = n x 1 vector, contains ones for nodes in region and zeros
%                 elsewhere 
% covOut        = p x p covariance matrix, where p is the number of nodes in the
%                 selected region

regionIndex=find(nodeIDoutRect);%% finds selected nodes
lengthRegion=length(regionIndex);

covOut=zeros(lengthRegion,lengthRegion);

[ind1,ind2]=meshgrid(regionIndex,regionIndex);                  % creates a grid with selected nodes as axis
[ind3,ind4]=meshgrid(1:lengthRegion,1:lengthRegion);            % creates a grid with length of new region as axis
ind1=reshape(ind1,[],1);
ind2=reshape(ind2,[],1);
ind3=reshape(ind3,[],1);
ind4=reshape(ind4,[],1);
covOut(sub2ind([lengthRegion,lengthRegion],ind3,ind4))...
    =covOrig(sub2ind(size(covOrig),ind1,ind2));                       % relates old marix V to new matrix tempV