% snippet to calculate in and out of tolerance percentage

close all;
clear ;
% dbstop if error;
run 'C:\Users\babu_m\Documents\GitHub\source\initFiles.m';
% path(path,'\\megara\dlmr\Manoj\all data sensitivity inner'); % all data files in megara

eps         =10E-16;
fem         =femInit();
fileName{1} =strcat(inpFiles,'\inner504_Coarse.inp');
fem         =importMultiMesh(fem, fileName);%has provisions for domain separation within code
fem         =femPreProcessing(fem);
domainID    =1;

coordi          =fem.xMesh.Node.Coordinate;
nodeIdDomain    =fem.Domain(domainID).Node;
nodeCoord       =fem.xMesh.Node.Coordinate(nodeIdDomain,:); % coordinates for the domain
%% load pre-defined part regions
partRegions     =load('inner9regions504CoarseMesh.mat');% contains variable nodeIDoutRect
nodeIDoutRectDense   =partRegions.nodeIDoutRect;

devPatterns		=load('simAutoCorDev504InnerBatches8NormDev.mat');%  
devPatterns		=devPatterns.simData;

seqPred			=load('504InnWokrig2batch45basis_1.0Tol_0.0W2.mat');%
seqPred			=seqPred.seqPred;

% loading gpr predictions
load('gprPredInner504_9regions7Batch.mat');

%% loading key points from predefined selection
keyPointMeshFile= 'selectedKeyInnerL504coarse.stl';
[Mnp,idMnp]		=getKeyPoints4mMesh(keyPointMeshFile,domainID,nodeCoord);
%% finding key points in each region
keyPointIndex			= zeros(size(nodeCoord,1),1);
keyPointIndex(idMnp)	= 1;
nodeIDoutRectCoarse		= bsxfun(@times,keyPointIndex,nodeIDoutRectDense);
%%
devTest     =devPatterns(7).FemDevDomain;
devKey		=devTest(:,idMnp);
nBasis		=45;%
maxSnap		=9;
numberCompleteMeasure =2;
Nparts		=size(devTest,1);
tol			=1;
cutoff      =.95; % prob above which tolerance predicition is taken as true value

for partNum = numberCompleteMeasure+1:Nparts
	
	
	for snapNo=1:maxSnap
		ytt		=seqPred(partNum).Snap(snapNo).DevKey;
		VarYtt	=seqPred(partNum).Snap(snapNo).VarKey;
		[allCorrect,~]	=findInOutTol(devKey(partNum,:)',ytt',VarYtt,tol,cutoff);
		perAllCorrect(partNum,snapNo)	=sum(allCorrect)/length(allCorrect);
		
		% for gpr
		yttSpace		=gprPred(partNum).Snap(snapNo).DevKey;
		VarYttSpac		=gprPred(partNum).Snap(snapNo).VarDevKey;
		[allCorrect,~]	=findInOutTol(devKey(partNum,:)',yttSpace',VarYttSpac,tol,cutoff);
		perAllCorrGPR(partNum,snapNo)	=sum(allCorrect)/length(allCorrect);
		
	end
end

perAllCorrect(1:numberCompleteMeasure,:)=[];
meanPercentST=mean(perAllCorrect);
perAllCorrGPR(1:numberCompleteMeasure,:)=[];
meanPercentS=mean(perAllCorrGPR);