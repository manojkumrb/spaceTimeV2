function [devRegion,varRegion]=getRegionDevVarInner(nodeIDoutRectCoarse,regionID,yt,VarYt,iMnp)
%
% Gets deviations and variance corresponding to the key oints in the region being measured
% from the deviation information of the keypoints for whole part
% remaps yt to global and then picks key points corresponding to a given
% region
%
% nodeIDoutRect     = no of nodes X no of regions matrix, with ones
%                     to represent the node belonging to a region
% yt                = deviation vector of all nodes in the part
% varYt             = variance vector of all nodes in the part
% devRegion         = deviation vector of all nodes in the region
% VarRegion         = variance vector of all nodes in the region


regionIndex=find(nodeIDoutRectCoarse(:,regionID));%% finds selected nodes
allDev=zeros(size(nodeIDoutRectCoarse,1),1);
allDev(iMnp)=yt;
allVar=zeros(size(nodeIDoutRectCoarse,1),1);
allVar(iMnp)=VarYt;
devRegion=allDev(regionIndex);
varRegion=allVar(regionIndex);

end