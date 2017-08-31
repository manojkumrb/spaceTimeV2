% get performance of combination of all possible deterministic sequences and saving the values for
% further evaluation.
close all;
clear ;
% dbstop if error;
run 'C:\Users\babu_m\Documents\GitHub\source\initFiles.m';
path(path,'\\megara\dlmr\Manoj\all data sensitivity inner'); % all data files in megara

eps         =10E-16;
fem         =femInit();
fileName{1} =strcat(inpFiles,'\door_inner.inp');%\top_hat.inp');%\halo.inp');%
% fem       =importMesh(fem,fileName{1});%'top_hat.inp'); %use only for top hat
fem         =importMultiMesh(fem, fileName);%has provisions for domain separation within code
fem         =femPreProcessing(fem);
domainID    =1;
idPart      =domainID;

coordi          =fem.xMesh.Node.Coordinate;
nodeIdDomain    =fem.Domain(idPart).Node;
nodeCoord       =fem.xMesh.Node.Coordinate(nodeIdDomain,:); % coordinates for the domain
%% load pre-defined part regions
partRegions     =load('inner9Regions.mat');% contains variable nodeIDoutRect
nodeIDoutRectDense   =partRegions.nodeIDoutRect;

devPatterns		=load('simAutoCorDevInnerBatchesCombined.mat');%  %hingeDevArSimLoc
devPatterns		=devPatterns.simData;
devPatterns		=devPatterns(8).FemDevDomain; % as the 6th batch is used for all tests

%% loading key points from predefined selection
selNodes=load('doorInnerSelNodes.mat');
selNodes=selNodes.selNodes;
Mnp=selNodes.Coord;
iMnp=selNodes.ID;

devCenterd			=devPatterns;%devCenterd1+randn(size(devPatterns))*10e-3;%devPatterns;x=z*sigma+mu;; % set equal to input in a coarse mesh
devCenterdKeypoint  =devCenterd(:,iMnp);

%% finding key points in each region
keyPointIndex		  = zeros(size(nodeIDoutRectDense,1),1);
keyPointIndex(iMnp)	  = 1;
nodeIDoutRectCoarse= bsxfun(@times,keyPointIndex,nodeIDoutRectDense);

%% load optimised error process covariance parameters
load('optimisedHypParmInner.mat'); % hypInnerFitc.mat');%contains struct hyp in gp function format(i.e. hyp.cov=log(parameters)) for matern covariance

sigmaMes    =1E-5;
nu.Type     ='diag';   % it turned out that, diagonal small scale covariance is better than a structured covariance function
nu.StdDev   =1e-4;
devAll      =devCenterd;
nBasis		=35;%20:10:60;    %     % set containing number of basis vectors

%% kalman recursion
	[keyEigVec,R,V,H,covErr]=getSystemMatricesSampled(nodeCoord,devAll,...
		sigmaMes,nu,iMnp,nBasis);
	devKey=devAll(:,iMnp);
	
	%% eigen interpolation
	fileNameCoarse= 'innerSelNodes - Copy.inp'; % the mesh of key points created fro cgal
	interpEigVec=getEigenInterp(nodeCoord,fileNameCoarse,keyEigVec, domainID);
	% interpEigVec=load('interpEigVecInner3_400.mat');
	% interpEigVec=interpEigVec.interpEigVec;
	%%
	zt          =devKey; % measurements
	tol         =2;    % is a tricky value affects the adaptivity a lot..currently works well when tol is set to 2
	countPart   =1;
	maxSnap     =9;
	numberCompleteMeasure= 10; % to stabalize the prediction
	pc          =devKey*keyEigVec;
	
	% at time zero
	a00         =pc(1,:)';
	p00         =eye(size(keyEigVec,2));
	ytp1        =zeros(size(devKey,1),size(Mnp,1));
	ytt         =zeros(size(ytp1));
	VarYtt      =zeros(size(ytp1));
	VarYtp1     =zeros(size(ytp1));
	
	%% initialization for dense mesh
	yttU        =zeros(size(devAll));
	ytp1U       =zeros(size(devAll));
	VarYttU		=zeros(size(devAll));

% to stabalise the filtering equations
for i=1:numberCompleteMeasure
	
		if i==1
			attm1=a00;
			pttm1=p00;
		else
			% predicting the state variables and their variances for (t+1)
			[attm1,pttm1]=getKalmanStatettm1(ptt,att,H,covErr);
		end
		
		%observation
		[kt,att,ptt]=getKalmanStatett(pttm1,attm1,keyEigVec,R,V,zt(i,:)');
		% prediction
		[ytt(i,:),VarYtt(i,:),ytp1(i,:),VarYtp1(i,:)]=getPredictions(keyEigVec,att,ptt,H,V,R,covErr);
	
end

%% deterministic pattern search

regionsMes.Seq(1).Pattern(1).Dev=[];

for nMes=3:9
	
	% creating all possible combinations as measuring regions 1,5,3 and 1,3,5
	% give the same result
	cTotal = nchoosek(1:9,nMes);

	
	for kk=1:size(cTotal,1)
		
		for i=numberCompleteMeasure+1:size(zt,1) % i for each part
			
			allMesReg=cTotal(kk,:); % all the measured regions
					
			
			[tempV,tempR,tempEigVec,tempZ,regionIndex]=getRegionSysMatInner(nodeIDoutRectCoarse,...
				allMesReg,V,R,interpEigVec,devAll(i,:),iMnp);
			
			% updating the state variables after partial measurement
			[tempKt,att,ptt]=getKalmanStatett(pttm1,attm1,tempEigVec,tempR,tempV,tempZ');
			
			% predicting deviations for the whole part using updated state variables
			[ytt(i,:),VarYtt(i,:),ytp1(i,:),VarYtp1(i,:),yttU(i,:),ytp1U(i,:),VarYttU(i,:)]=getPredictionsNewInner(interpEigVec,...
				keyEigVec,tempEigVec,att,ptt,H,V,R,covErr,hyp,nodeCoord,...
				nodeCoord(iMnp,:),nodeCoord(regionIndex,:),tempZ);
			

% 			RMSE calculation for the step
			diff=bsxfun(@minus,devCenterdKeypoint(i,:),ytt(i,:));
			rmseT=sqrt(sum(diff.^2)/size(devCenterdKeypoint,2));  % rmse for key points
% 			diff=bsxfun(@minus,devAll(i,:),yttU(i,:));
% 			rmseT(nSnaps)=sqrt(sum((diff).^2)/size(devAll,2));  % rmse for all points
			
			
			% saving predictions
			regionsMes.Seq(kk).Regions=allMesReg;
			regionsMes.Seq(kk).Pattern(i).Dev=yttU(i,:);
			regionsMes.Seq(kk).Pattern(i).Var=VarYttU(i,:);
			regionsMes.Seq(kk).Pattern(i).Rmse=rmseT;
		
		end
		fprintf('%d of %d sequences completed with %d regions \n',kk,size(cTotal,1),nMes);
		
		
	end
	fileName=sprintf('detAllSeq%dReg',nMes);
	save(fileName,'regionsMes','-v7.3')

	
end
