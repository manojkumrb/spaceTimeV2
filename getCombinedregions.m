function allPoints=getCombinedregions(nodeIDoutRect,allMesReg)
% nodeIDoutRect             = no of nodes X no of regions matrix, with ones
%                             to represent the node belonging to a region
% allMesReg                 = list of all regions measured, numbers
%                             correspond to columns of nodeIDoutRect

allPoints=zeros(size(nodeIDoutRect,1),1);

for j=1:length(allMesReg)
    allPoints=allPoints|nodeIDoutRect(:,allMesReg(j));
end

end