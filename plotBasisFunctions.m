% a separate code to create eigen basis representaion for illustrations

close all;
clear ;
dbstop if error;

run 'C:\Users\babu_m\Documents\GitHub\source\initFiles.m';
path(path,'\\megara\dlmr\Manoj\all data sensitivity inner'); % all data files in megara

eps         =10E-16;
fem         =femInit();
fileName{1} =strcat(inpFiles,'\door_inner.inp');%\top_hat.inp');%\halo.inp');%
% fem       =importMesh(fem,fileName{1});%'top_hat.inp'); %use only for top hat
fem         =importMultiMesh(fem, fileName);%has provisions for domain separation within code


%% SIMULATED autocorrelated part variation
% inner
devPatterns		=load('simAutoCorDevInnerBatchesCombined.mat');%  %hingeDevArSimLoc
devPatterns		=devPatterns.simData;
devTrain			=devPatterns(8).FemDevDomain;
devTest				=devPatterns(8).FemDevDomain;
% loading key points from predefined selection
selNodes=load('doorInnerSelNodes.mat');
selNodes=selNodes.selNodes;
Mnp=selNodes.Coord;
iMnp=selNodes.ID;
% flatplate
load('flatPlateDevArLocators.mat'); % contains variable femDevArFlatPlate



fem         =femPreProcessing(fem);
domainID    =1;

coordi          =fem.xMesh.Node.Coordinate;
nodeIdDomain    =fem.Domain(domainID).Node;
nodeCoord       =fem.xMesh.Node.Coordinate; % coordinates for the domain




% some inputs necessary for the function
nBasis		= 12;
sigmaMes    =1E-4;
nu.Type     ='diag';
nu.StdDev   =1e-4;

keyNodeIds	=nodeIdDomain; % in case of large parts


%% for hinge
% devAll			=zeros(size(femDevArFlatPlate,1),size(coordi,1));
% devAll(:,nodeIdDomain)= femDevArFlatPlate;%dev;%
% 
% [eigVec,R,V,H,covErr]=getSystemMatricesSampled(nodeCoord,devAll,...
% 	sigmaMes,nu,keyNodeIds,nBasis);

% %  Plotting eigen basis
% for i=1:size(eigVec,2)
% 	contourDomainPlot(fem,domainID,eigVec(:,i),0);
% 	export_fig(sprintf('pattern%i',i),'-png','-r400','-transparent');
% 	contourDomainPlot(fem,domainID,femDevArFlatPlate(i,:),0);
% 	export_fig(sprintf('Basis%i',i),'-pdf','-transparent');
% 	export_fig(sprintf('Basis%i',i),'-png','-r400','-transparent');
% 		print(sprintf('Basis%i',i),'-dpng')
% end


%% for inner

[keyEigVec,R,V,H,covErr]=getSystemMatricesSampled(nodeCoord,devTrain,...
sigmaMes,nu,iMnp,nBasis);
% eigen interpolation
fileNameCoarse= 'innerSelNodes - Copy.inp';
interpEigVec=getEigenInterp(nodeCoord,fileNameCoarse,keyEigVec, domainID);

%  Plotting eigen basis
for i=1:2%size(eigVec,2)
	contourDomainPlot(fem,domainID,interpEigVec(:,i),0);
	export_fig(sprintf('pattern%i',i),'-png','-r400','-transparent');
	contourDomainPlot(fem,domainID,devTrain(i,:),0);
% 	export_fig(sprintf('Basis%i',i),'-pdf','-transparent');
	export_fig(sprintf('Basis%i',i),'-png','-r400','-transparent');
	% 	print(sprintf('Basis%i',i),'-dpng')
end

