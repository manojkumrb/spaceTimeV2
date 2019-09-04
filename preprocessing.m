% code for all the prepocessing to the mesh files
% select regions,% select key points
% for key point sensitivity check keyPointSensitivity.m

close all; clear ;

run 'C:\Users\babu_m\Documents\GitHub\source\initFiles.m';

fem=femInit();
fileName{1}='C:\Users\babu_m\Documents\GitHub\source\input\inner504_Coarse.inp';%

fem=importMultiMesh(fem, fileName);
fem=femPreProcessing(fem);

domainID=1;
nodeIdDomain=fem.Domain(domainID).Node;
nodeCoord=fem.xMesh.Node.Coordinate(nodeIdDomain,:); % coordinates for the domain
nRegions=9;


%% selecting regions
% nodeIDoutRect = getRegionsGui(fem,domainID,nRegions);
% save('inner9regions405CoarseMesh.mat','nodeIDoutRect');

%% plotting point cloud on the mesh
% cloud = load('C:\Users\babu_m\Desktop\output_VAM_05_09_2017_08_18_22_0.vtx');
% % plotCloudOnMesh(fem,cloud,2,0.1)


%% graphical selction of key points
[ ~, ~,iMnp  ] = uniformKeyPointSelection( fem,domainID,20);
[boundNodeID]=getBoundaryNodes(fem,domainID, 150); % doesnot work as intended for new inner at the flage
% has very small elements

% picking the edge node graphically to key points in areas not automatically selsected
% meshplotAxisDefined(fem,domainID);
% fig=gcf;
% title('Select edge node');
% ax=gca;
% ax.Clipping='off'; % prevent clipping when zooming matlab
% [nodeIDsel,~]=guiSelection(fig,ax,nodeCoord);

% manualSelect=[10542;10546;10549;10552;10554;10557;22201;2274;22199;...
% 	22194;22191;21865;21867;21870;21873;21876;21879;21881;21884;21886;21889;...
% 	20821;20825;20828;20832;20836;21029;20126;21022;];
% iMnp=[iMnp;manualSelect];

iMnp=[iMnp;boundNodeID];
uniqeMnp= unique(iMnp);
Mnp=nodeCoord(uniqeMnp,:);
selNodes.ID=uniqeMnp;
selNodes.Coord=Mnp;

% save('keyPt405InnerCoarse20.mat','selNodes');


%% loading from predefined selection
selNodes=load('keyPt405InnerCoarse15.mat');
selNodes=selNodes.selNodes;

Mnp=selNodes.Coord;
iMnp=selNodes.ID;
% iMnp=dsearchn(nodeCoord,Mnp);
%% checking selected nodes
meshplotAxisDefined(fem,1);
hold all
scatter3(Mnp(:,1),Mnp(:,2),Mnp(:,3),'*');

%% opening stl to copy node id and coordinates



%%
devPatterns= load('devAutoCorrInner_dense.mat');%flatPlateDevArLocators.mat');%
devPatterns=devPatterns.femDevDomain;%femDevArFlatPlate;
devCenterd1=bsxfun(@minus,devPatterns,mean(devPatterns)); % to calculate the rmse, used nowhere else
devSelected= devPatterns(:,iMnp);
devCenterd=bsxfun(@minus,devSelected,mean(devSelected));

V=load('eigVecInner3_400.mat');
V=V.V;
pc=devCenterd*V;%plot(pc);
x=pc*V';
diff=devCenterd-x;
figure
imagesc(diff);
