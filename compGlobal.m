%% a file to including loops for batches, basis, weights
% removed all redundat code, refer compareEigInnerNewAdapCriter for code to
% plot and other redundancies

% {
% The block comment line to enable effective debug -remove the space after '%'
%close all;
clear ;
dbstop if error;

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
%}
%% load pre-defined part regions
partRegions     =load('inner9Regions.mat');% contains variable nodeIDoutRect
nodeIDoutRectDense   =partRegions.nodeIDoutRect;

%% SIMULATED autocorrelated part variation
devPatterns		=load('simAutoCorDevInnerBatchesCombined.mat');%  %hingeDevArSimLoc
devPatterns		=devPatterns.simData;

%% loading key points from predefined selection
selNodes=load('doorInnerSelNodes.mat');
selNodes=selNodes.selNodes;
Mnp=selNodes.Coord;
iMnp=selNodes.ID;

%% finding key points in each region
keyPointIndex		  = zeros(size(nodeIDoutRectDense,1),1);
keyPointIndex(iMnp)	  = 1;
nodeIDoutRectCoarse= bsxfun(@times,keyPointIndex,nodeIDoutRectDense);

%% load optimised error process covariance parameters
load('hypInnerFitc.mat');% optimisedHypParmInner.mat'); %contains struct hyp in gp function format(i.e. hyp.cov=log(parameters)) for matern covariance
sigmaMes    =1E-4;
nu.Type     ='diag';
nu.StdDev   =1e-4;

%% setting test and train data and other parameters
testBatch		=8;
nTrainBatch		=[1,2,3,4,5,6];
nBasis			=35;%5:5:70;    %     % set containing number of basis vectors
infoWeight		= .8;
distweight		= .2;

devTest				=devPatterns(testBatch).FemDevDomain;
devKey  =devTest(:,iMnp);


for j=1:length(nTrainBatch)
	
	
	devTrain			=[];
	
	for l=1:nTrainBatch(j)        % adding train data
		devTrain		=[devTrain;devPatterns(l).FemDevDomain];
	end
	
	
	%%
	for w=1:length(infoWeight)
		
		for ww=1:length(distweight)
			
			for k=1:length(nBasis)
				
				%% kalman recursion
				[keyEigVec,R,V,H,covErr]=getSystemMatricesSampled(nodeCoord,devTrain,...
					sigmaMes,nu,iMnp,nBasis(k));
				
				
				
				%% eigen interpolation
				fileNameCoarse= 'innerSelNodes - Copy.inp';
				interpEigVec=getEigenInterp(nodeCoord,fileNameCoarse,keyEigVec, domainID);
				% interpEigVec=load('interpEigVecInner3_400.mat');
				% interpEigVec=interpEigVec.interpEigVec;
				%%

				tol         =4;    % is a tricky value affects the adaptivity a lot..currently works well when tol is set to 2
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
				yttU        =zeros(size(devTest));
				ytp1U       =zeros(size(devTest));
				VarYttU		=zeros(size(devTest));
				
				%% to stabalise the filtering equations
				for i=1:numberCompleteMeasure
					
					if i==1
						attm1=a00;
						pttm1=p00;
					else
						% predicting the state variables and their variances for (t+1)
						[attm1,pttm1]=getKalmanStatettm1(ptt,att,H,covErr);
					end
					
					%observation
					[kt,att,ptt]=getKalmanStatett(pttm1,attm1,keyEigVec,R,V,devKey(i,:)');
					% prediction
					[ytt(i,:),VarYtt(i,:),ytp1(i,:),VarYtp1(i,:)]=getPredictions(keyEigVec,att,ptt,H,V,R,covErr);
					
				end
				
				%% Adaptive and partial measurement
				tic
				seqPred(size(devKey,1)).Snap(1).Dev=[];
				rmseTp1=zeros(size(devKey,1),1);
				
				
				for i=numberCompleteMeasure+1:size(devKey,1) % i for each part
					
					nSnaps=1;
					allMesReg=[];
					nRegions=1:size(nodeIDoutRectCoarse,2);
					rmseT=zeros(maxSnap,1);
					selRegion=zeros(maxSnap,1);
					
					while ~isempty(nRegions)&&(nSnaps<=maxSnap)
						selIndex=zeros(length(nRegions),1);
						distance=zeros(length(nRegions),1);
						
						% prediction based on t-1 to predict first region to measure
						for jj=1:length(nRegions)
							if nSnaps==1
								% when no region is measured use previous prediction
								% (ytp1) and distance is set to zero
								devPart=ytp1(i-1,:);
								varPart=VarYtp1(i-1,:);
								distance(jj)=0;
							else
								% when some regions are measured use updated prediction (ytt)
								devPart=ytt(i,:);
								varPart=VarYtt(i,:);
								[distance(jj)]=getDistBwRegions(nodeIDoutRectCoarse,nodeCoord, nRegions(jj),selRegion(nSnaps-1));
							end
							
							[devRegion,varRegion]=getRegionDevVarInner(nodeIDoutRectCoarse,nRegions(jj),devPart,varPart,iMnp);
							% calculating the probability of being out of tolerance
							[selIndex(jj)]=getRegionScore(tol, devRegion, varRegion);
							
						end
						
						combinedScore=getCombinedScore(selIndex,distance,infoWeight(w),distweight(ww));
						
						[~,indexMax]=max(combinedScore);
						
						selRegion(nSnaps)=nRegions(indexMax);
						allMesReg=[allMesReg;selRegion(nSnaps)]; % all the measured regions
						
						[tempV,tempR,tempEigVec,tempZ,regionIndex]=getRegionSysMatInner(nodeIDoutRectCoarse,...
							allMesReg,V,R,interpEigVec,devTest(i,:),iMnp);
						
						% updating the state variables after partial measurement
						[tempKt,att,ptt]=getKalmanStatett(pttm1,attm1,tempEigVec,tempR,tempV,tempZ');
						
						% predicting deviations for the whole part using updated state variables
						[ytt(i,:),VarYtt(i,:),ytp1(i,:),VarYtp1(i,:),yttU(i,:),ytp1U(i,:),VarYttU(i,:)]=getPredictionsNewInner(interpEigVec,...
							keyEigVec,tempEigVec,att,ptt,H,V,R,covErr,hyp,nodeCoord,...
							nodeCoord(iMnp,:),nodeCoord(regionIndex,:),tempZ);
						% 			contourDomainPlot(fem,1,mpredictedU,1);
						
						
						% removing the measured region from the list
						nRegions(indexMax)=[];
						% 			RMSE calculation for the step
						diff=bsxfun(@minus,devKey(i,:),ytt(i,:));
						rmseT(nSnaps)=sqrt(sum(diff.^2)/size(devKey,2));  % rmse for key points
						% 			diff=bsxfun(@minus,devAll(i,:),yttU(i,:));
						% 			rmseT(nSnaps)=sqrt(sum((diff).^2)/size(devAll,2));  % rmse for all points
						
						% saving sequential predictions
						seqPred(i).Snap(nSnaps).Dev=yttU(i,:);
						seqPred(i).Snap(nSnaps).Var=VarYttU(i,:);
						seqPred(i).Snap(nSnaps).RmseT=rmseT;
						seqPred(i).Snap(nSnaps).Region=allMesReg;%selRegion;
						nSnaps=nSnaps+1;
					end
					
					% RMSE for (t+1) prediction
					diff=bsxfun(@minus,devKey(i,:),ytp1(i-1,:));
					rmseTp1(i)=sqrt(sum((diff).^2)/size(devKey,2));  % rmse without noisy input
					seqPred(i).Snap(maxSnap).rmseTp1=rmseTp1(i);
					seqPred(i).Snap(maxSnap).Ytp0=ytp1(i-1,:);
					
					% predicting the state variables and their variances for (t+1)
					[attm1,pttm1]=getKalmanStatettm1(ptt,att,H,covErr);
					
					
				end
				
				toc
				%% saving data
				fileString=sprintf('innerEigWithkrigYtp1%ibatch%ibasis_%.1fW1_%.1fW2.mat',nTrainBatch(j),nBasis(k),infoWeight(w),distweight(ww));
				save(fileString,'seqPred','-v7.3');
				clear seqPred;
				
			end
			
			
		end
		
		
	end
	
	
end
%%
% plotting rmse sensitivity

% nBasis		=40;%[5:5:60];
numberCompleteMeasure=10;
maxSnap=9;

% plot parameters
markers = {'o','s','d','^','v','x','+','*','.','>','<','p','h','o','s','d','^','v','x'};
cmapRmsePlot = cbrewer('qual','Set1',length(nBasis));

% reading data from file

for j=1:length(nTrainBatch)
	
	for w=1:length(infoWeight)
		
		for ww=1:length(distweight)
			
			for k=1:length(nBasis)
				
				fileString=sprintf('innerEigWithkrigYtp1%ibatch%ibasis_%.1fW1_%.1fW2.mat',nTrainBatch(j),nBasis(k),infoWeight(w),distweight(ww));
				seqPred=load (fileString);
				seqPred=seqPred.seqPred;
				
				for i=numberCompleteMeasure+1:length(seqPred)
					rmsePlot(:,i)=seqPred(i).Snap(maxSnap).RmseT;		% accessing each of the 50 replication
				end
				
				rmseAvg(:,k)=sum(rmsePlot,2)./(length(seqPred)-(numberCompleteMeasure+1)); % averaging all replications
				
			end
			
			% plotting eig cmparison
			lineRmseEig=plot(rmseAvg);
			hold on;
			
		end
		
		
	end
	
	
end

load('S:\DLMR\Manoj\all data sensitivity inner\bestDeterministic.mat')
plot(rmseAvgDet)



