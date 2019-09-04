%% usage contourDomainPlot(fem,domainID,userData,colourBar,ax)
% creating contourplot for a given domain within fem structure, for data
% specified by the user, length(userdata)== number of nodes in the domain.
% colourBar = 0 or 1, creates a colourBar when = 1
% if ax is empty a new axix is created

function contourDomainPlot(fem,domainID,userData,colourBar,varargin)

noOfNodesDomain=length(fem.Domain(domainID).Node);
if length(userData)~= noOfNodesDomain
    display('user data doesnot match the number of nodes in domain');

    return
end

nodeIDs=fem.Domain(domainID).Node;
if isempty(varargin)
    figure;
    ax=gca;
else
    ax=varargin{1};
end
% title('Predicted Deviation','fontweight','bold','fontsize',13);
% for all domains together for contourplot
AllNodesWholeDomain=zeros(size(fem.xMesh.Node.Coordinate,1),1);
AllNodesWholeDomain(nodeIDs)=userData;

% hold all
fem.Post.Options.ParentAxes=ax;
fem.Sol.UserExp=AllNodesWholeDomain';
fem.Post.Contour.Domain=domainID;
fem.Post.Contour.ContourVariable='user';
fem.Post.Contour.Resolution=1;
fem.Post.Contour.MinRangeCrop =-inf;
fem.Post.Contour.MaxRangeCrop=inf;
contourPlot(fem)

if colourBar==1
    barHandle=colorbar(ax);
    barHandle.FontSize=12;
%     set(barHandle, 'Position', [.95 .11 .025 .75])
    set(get(barHandle,'Ylabel'),'string',...
        'surface normal deviation (mm)','fontweight','bold','fontsize',13)    
else
     barHandle=colorbar;
     delete(barHandle);
end
axes(ax)
axis equal
axis tight
ax.Clipping = 'off'; 
axis off
view([0 0]);
