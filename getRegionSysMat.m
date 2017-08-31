function [tempV,tempR,tempEigVec,tempZ,regionIndex]=getRegionSysMat(nodeIDoutRect,allMesReg,V,R,eigVec,zt)

%
% nodeIDoutRect             = no of nodes X no of regions matrix, with ones
%                             to represent the node belonging to a region
% allMesReg                 = list of all regions measured, numbers
%                             correspond to columns of nodeIDoutRect
% eigVec                    = eigen basis (EOF) for all nodes 
% zt                        = Deviation information for all nodes
% tempV,tempR,tempEigVec,
% tempZ                     = V,R,eigVec and Z corresponding to the
%                             selected regions
% regionIndex               = domain node IDs of all selected nodes

allPoints=zeros(size(nodeIDoutRect,1),1);

for j=1:length(allMesReg)
    allPoints=allPoints|nodeIDoutRect(:,allMesReg(j));
end

tempV=rearrangeCovarianceMat(V,allPoints);
tempR=rearrangeCovarianceMat(R,allPoints);
regionIndex=find(allPoints);%% finds selected nodes
tempEigVec=eigVec(regionIndex,:);
tempZ= zt(regionIndex);
end 