% calculate and plot deviations from Nominal using COP with options for
% plotting as subplots and also cop and normal overlay with mesh
% {
close all;
clear ;
dbstop if error;

init		= 'C:\Users\babu_m\Documents\GitHub\source\initFiles.m';
run(init); 
% path(path,'\\megara\dlmr\Manoj\all data sensitivity inner'); % all data files in megara

eps         =10E-16;
fem         =femInit();
fileName{1} =strcat(inpFiles,'\inner405_Coarse.inp');
fem         =importMultiMesh(fem, fileName);%has provisions for domain separation within code
fem         =femPreProcessing(fem);
domainID    =1;

coordi          =fem.xMesh.Node.Coordinate;
nodeIdDomain    =fem.Domain(domainID).Node;
nodeCoord       =fem.xMesh.Node.Coordinate(nodeIdDomain,:); % coordinates for the domain
normal			=fem.xMesh.Node.Normal(nodeIdDomain,:);

%}

load ('visNodeIDoutRect.mat');
%%
count=0;
devMeshNodeCombined		=zeros(size(nodeCoord,1),1);
for i=1:9
    count=count+1;
	file='C:\Users\babu_m\Documents\GitHub\vnv\ExportedTiles\Door_inner_single_panel_CoP_marked_SUBSAMPLED.vtx';
%     file=sprintf('C:\\Users\\babu_m\\Documents\\GitHub\\vnv\\ExportedTiles\\%d.xyz',i);%Halo_part %d.txt',i)%Hinge Front Right_part %d.txt',i);% 
    cloud=load(file); % contains variable dev
%     devMeshNode(:,count)=getNormalDevPoints2Points(nodeCoord,normal,cloud,12,1);    % normal search dist, area of foot of cylinder
    devMeshNodeAll=getNormalDevPoints2Points(nodeCoord,normal,cloud,12,1);    % normal search dist, area of foot of cylinder

	uniqueID	= find(devMeshNode(:,count));
	devMeshNodeCombined(uniqueID)	= devMeshNode(uniqueID,count);

end



% to check the overlap of cop and mesh files
percentage=0.1;
nodeNormalSample(fem,domainID,percentage)
plotCloudOnMesh(fem,cloud,domainID,percentage)
%
figure
for i=count%:size(devHatLower,2)
    % title('Mean Deviation','fontweight','bold')
    ax=axes;%
    hold all
    fem.Post.Options.ParentAxes=ax;
    fem.Post.Contour.Domain=domainID;
    fem.Post.Contour.ContourVariable='user';
    allNodesDev=zeros(size(fem.xMesh.Node.Coordinate,1),1);
    allNodesDev(nodeIDs)=devMeshNode(:,i);
    fem.Sol.UserExp=allNodesDev';
    fem.Post.Contour.Resolution=1;
    contourPlot(fem)
    barHandle=colorbar;
    set(get(barHandle,'Ylabel'),'string','Deviations in mm','fontweight','bold')
    axis off
end

toc