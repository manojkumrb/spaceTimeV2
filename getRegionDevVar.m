function [devRegion,varRegion]=getRegionDevVar(nodeIDoutRect,regionID,yt,VarYt)
%
% Gets deviations and variance corresponding to the region being measured
% from the deviation information of the whole part
%
% nodeIDoutRect     = no of nodes X no of regions matrix, with ones
%                     to represent the node belonging to a region
% yt                = deviation vector of all nodes in the part
% varYt             = variance vector of all nodes in the part
% devRegion         = deviation vector of all nodes in the region
% VarRegion         = variance vector of all nodes in the region


regionIndex=find(nodeIDoutRect(:,regionID));%% finds selected nodes
devRegion=yt(regionIndex);
varRegion=VarYt(regionIndex);

end