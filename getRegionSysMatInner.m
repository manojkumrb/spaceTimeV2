function [tempV,tempR,tempEigVec,tempZ,regionIndex,keyPtLocIndex]=getRegionSysMatInner(nodeIDoutRectCoarse,allMesReg,V,R,eigVec,zt,iMnp)

%
% nodeIDoutRectcoarse       = no of nodes X no of regions matrix, with ones
%                             to represent the key points belonging to a region,
%                             in global coordinates, ie length = length of
%                             total nodes in the domain (source/dense mesh)
% allMesReg                 = list of all regions measured, numbers
%                             correspond to columns of nodeIDoutRect
% eigVec                    = eigen basis (EOF) for all nodes 
% zt                        = Deviation information for all nodes
% iMnp						= mesh node id of key points in global 
% tempV,tempR,tempEigVec,
% tempZ                     = V,R,eigVec and Z corresponding to the
%                             selected regions
% regionIndex               = domain node IDs of all selected nodes


nodeIDoutRectCoarseLocal	= nodeIDoutRectCoarse(iMnp,:); % converts to local indexing (counts from 1: length(iMnp))
allPointsLocal				=zeros(size(nodeIDoutRectCoarseLocal,1),1);
allPoints					=zeros(size(nodeIDoutRectCoarse,1),1);

for j=1:length(allMesReg)
	allPoints		=allPoints|nodeIDoutRectCoarse(:,allMesReg(j));
    allPointsLocal	=allPointsLocal|nodeIDoutRectCoarseLocal(:,allMesReg(j));
end


% all points is the index of the points to be selected, should be in 
% local coordiantes, ie coordinates of the key points (too lazy to change the below function)
tempV	=rearrangeCovarianceMat(V,allPointsLocal);
tempR	=rearrangeCovarianceMat(R,allPointsLocal);


regionIndex		=find(allPoints);%% finds selected nodes
keyPtLocIndex	=find(allPointsLocal);
tempEigVec		=eigVec(regionIndex,:);
tempZ			=zt(regionIndex);
end 