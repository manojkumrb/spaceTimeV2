function [distance]=getDistBwRegions(nodeIDoutRectCoarse,nodeCoord, regionIdT,regionIdTm1)
% calcultes the distance between two region CGs

regionIndex		=find(nodeIDoutRectCoarse(:,regionIdT));%% finds selected nodes
keyPointsRegion    =nodeCoord(regionIndex,:);
regionIndexTm1		=find(nodeIDoutRectCoarse(:,regionIdTm1));%% finds selected nodes
keyPointsRegionTm1    =nodeCoord(regionIndexTm1,:);

cgT				= sum(keyPointsRegion)./size(keyPointsRegion,1);
cgTm1				= sum(keyPointsRegionTm1)./size(keyPointsRegionTm1,1);
distance	= sqrt(sum((cgT-cgTm1).^2));
end