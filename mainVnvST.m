%%  Code for running STAS methodology for the paper in JMS
%% Copyright 2019 <Authors>
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the "Software"), 
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be 
% included in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
% ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
% DEALINGS IN THE SOFTWARE.
%% This is a code which includes all the necessites to run the methodology as described in paper
% sensitivity analysis have to be run by modifying this code and spatial
% comparison needs to be replicated using GPML toolbox and the
% information in the paper


close all;
clear ;

load('femStrut.mat'); % loads the mesh file in a given format
% this is created using the VRM sofware from DLM group at WMG, University of
% Warwick
coordi          =fem.xMesh.Node.Coordinate;
nodeIdDomain    =fem.Domain(domainID).Node;
nodeCoord       =fem.xMesh.Node.Coordinate(nodeIdDomain,:); % coordinates for the domain

%% load pre-defined part regions
partRegions			=load('inner9regions504CoarseMesh.mat');% contains variable nodeIDoutRect nNodes X regions logical matrix
nodeIDoutRectDense  =partRegions.nodeIDoutRect;

%% SIMULATED autocorrelated part variation
load('simAutoCorDev504InnerBatches8NormDev.mat'); % simulated using VRM sftware
%% loading key points from predefined selection
load('keyPt405InnerCoarse20.mat');
Mnp=keyPoints.Coordiante;
idMnp=keyPoints.ID;
%% finding key points in each region
keyPointIndex			= zeros(size(nodeCoord,1),1);
keyPointIndex(idMnp)	= 1;
nodeIDoutRectCoarse		= bsxfun(@times,keyPointIndex,nodeIDoutRectDense);

%% load optimised error process covariance parameters
sigmaMes    =1E-12; % originally E-12

% coloured noise requires gp toolbox from Rasmussen
load('inner504HypOpt.mat');% %contains struct hyp in gp function format(i.e. hyp.cov=log(parameters)) for matern covariance
% hyp.lik		=log(0.00001);
% nu.Type     ='covFunc';%
% nu.CovPar   =hyp;

% error reduces when noise process is modelled as independent
% erors for simulation
nu.Type     ='diag';%
nu.StdDev   =1e-6;

%% setting test and train data and other parameters
nBasis		=45;% set containing number of basis vectors
tol         =1;    % design tolerance value for surface inspection
cutoff		=0.95;
countPart   =1;
maxSnap     =9;
numberCompleteMeasure= 2;

%%%%% test data input%%%%%%%%%%
devTest		= [simData(7).FemDevDomain];% matrix (repetitions X nodes)% 
devTrain	=[simData(1).FemDevDomain; simData(2).FemDevDomain];%simData(2).FemDevDomain ;% devSmooth';% 

%%%%%%%%%%%%%%%%%%%%%%%
devKey  =devTest(:,idMnp);
%% get system matrices
[~,R,V,H,covErr]=getSystemMatricesSampled(nodeCoord,devTrain,...
	sigmaMes,nu,idMnp,nBasis);

%% eigen interpolation
% interpEigVec	=getEigenInterp(nodeCoord, keyPointMeshFile,keyEigVec, domainID);
[~ ,~,eigVec]	=svds(devTrain,nBasis);
keyEigVec		=eigVec(idMnp,:);
interpEigVec	=eigVec;

%%
% to stabalize the prediction
pc          =devKey*keyEigVec;
nParts		=size(pc,1);
pcCent		=bsxfun(@minus,pc,mean(pc));
empStateCov	=(pcCent'*pcCent)./nParts;

empLoss		= devKey-(keyEigVec*pc')';
% at time zero
a00         =pc(1,:)';
p00         =1e10*eye(size(keyEigVec,2));%empStCov;%
ytp1        =zeros(1,size(Mnp,1));
ytt         =zeros(size(ytp1));
VarYtt      =zeros(size(ytp1));
VarYtp1     =zeros(size(ytp1));

%% initialization for dense mesh
yttU        =zeros(1,size(devTest,2));
ytp1U       =zeros(1,size(devTest,2));
VarYttU		=zeros(1,size(devTest,2));

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
	[ytt,VarYtt,ytp1,VarYtp1]=getPredictions(keyEigVec,att,ptt,H,V,R,covErr);	
end

%% Adaptive and partial measurement
tic
seqPred(size(devKey,1)).Snap(1).Dev=[];
rmseTp1=zeros(size(devKey,1),1);
countNonZer=1;
Nparts =20;%size(devTest,1);
for partNum		= numberCompleteMeasure+1:Nparts
	nSnaps=1;
	allMesReg=[];
	nRegions=1:size(nodeIDoutRectCoarse,2);
	rmseT=zeros(maxSnap,1);
	perAllCorrect=zeros(maxSnap,1);
	errorTStdDev=zeros(maxSnap,1);
	selRegion=zeros(maxSnap,1);	
	while ~isempty(nRegions)&&(nSnaps<=maxSnap)
		selIndex=zeros(length(nRegions),1);
		distance=zeros(length(nRegions),1);		
		% prediction based on t-1 to predict first region to measure
		for jj=1:length(nRegions)
			if nSnaps==1
				% when no region is measured use previous prediction
				% (ytp1) and distance is set to zero
				devPart=ytp1;
				varPart=VarYtp1;
				distance(jj)=0;
			else
				% when some regions are measured use updated prediction (ytt)
				devPart=ytt;
				varPart=VarYtt;
				[distance(jj)]=getDistBwRegions(nodeIDoutRectCoarse,nodeCoord, nRegions(jj),selRegion(nSnaps-1));
			end			
			[devRegion,varRegion]=getRegionDevVarInner(nodeIDoutRectCoarse,nRegions(jj),devPart,varPart,idMnp);
			% calculating the probability of being out of tolerance
			[selIndex(jj),regionData(jj).Entropy]=getRegionScore(tol, devRegion, varRegion);			
		end		
		regionIdBeforeSel =nRegions;		
		combinedScore=getCombinedScore(selIndex,distance,1,0);		
		[~,indexMax]=max(combinedScore);		
		selRegion(nSnaps)=nRegions(indexMax);
		allMesReg=[allMesReg;selRegion(nSnaps)]; % all the measured regions
		[tempV,tempR,tempEigVec,tempZ,regionIndex,keyPtLocIndex]=getRegionSysMatInner(nodeIDoutRectCoarse,...
			allMesReg,V,R,interpEigVec,devTest(partNum,:),idMnp);
		% updating the state variables after partial measurement
		[tempKt,att,ptt]=getKalmanStatett(pttm1,attm1,tempEigVec,tempR,tempV,tempZ');
		
		
		%% predicting deviations for the whole part using updated state variables
		[ytt,VarYtt,ytp1,VarYtp1,yttU,ytp1U,...
			VarYttU,VarYtp1Full]=getPredictionsNewInnerWOKrig(interpEigVec,...
			keyEigVec,tempEigVec,att,ptt,H,V,R,covErr,hyp,nodeCoord,...
			nodeCoord(idMnp,:),nodeCoord(regionIndex,:),tempZ);
	
		% setting measured points dev and variance
		ytt(keyPtLocIndex)		=devKey(partNum,keyPtLocIndex);
		VarYtt(keyPtLocIndex)	=0;
		
		

		%% removing the measured region  and clearing varibles from the list
		nRegions(indexMax)=[];
		
		%RMSE calculation for the step
		diff				=bsxfun(@minus,devKey(partNum,:),ytt);
		madT				=mean(abs(diff));
		mdT					=mean((diff));
		rmseT(nSnaps)		=sqrt(sum(diff.^2)/size(devKey,2));  % rmse for key points
		errorTVarDev(nSnaps)=var(diff);
		RHO2	= corr(devKey(partNum,:)',ytt','type','Spearman');
		bins	=100;
		
		[allCorrect,pYttSave]=findInOutTol(devKey(partNum,:)',ytt',VarYtt,tol,cutoff);
		perAllCorrect(nSnaps)	=sum(allCorrect)/length(allCorrect);

		% saving sequential predictions
			seqPred(partNum).Snap(nSnaps).Dev		=yttU;
			seqPred(partNum).Snap(nSnaps).DevKey	=ytt;
			seqPred(partNum).Snap(nSnaps).Var		=VarYttU;
			seqPred(partNum).Snap(nSnaps).VarKey	=VarYtt;
			seqPred(partNum).Snap(nSnaps).RmseT		=rmseT;
			seqPred(partNum).Snap(nSnaps).Region	=allMesReg;%selRegion;
			seqPred(partNum).Snap(nSnaps).RegionData=regionData;
			seqPred(partNum).Snap(nSnaps).RegionIdBeforeSel=regionIdBeforeSel;%selRegion;
			seqPred(partNum).Snap(nSnaps).MadT		=madT;
			seqPred(partNum).Snap(nSnaps).MdT		=mdT;
			seqPred(partNum).Snap(nSnaps).ErrorTKey	=diff;
			seqPred(partNum).Snap(nSnaps).Scor		=RHO2;
			seqPred(partNum).Snap(maxSnap).ErrorTVarDev	=errorTVarDev;
			

		%% plotting current predictions and errors for chosen given part
		% {
		if partNum==3
			
			if nSnaps ==1
				figComb=plotRmseErrorBarAndWrong(rmseT(1:nSnaps),sqrt(errorTVarDev(1:nSnaps)),perAllCorrect(1:nSnaps));
				
				panhandle	= uipanel('parent',figComb,'Position', [0 0.5 1 .5]);
				axOrigDev	= subplot(1,3,1,'Parent', panhandle,...
					'tag','axOrigDev', 'NextPlot', 'add');
				contourDomainPlot(fem,domainID,devTest(partNum,:),1,axOrigDev);
				origBarHandle	=colorbar(axOrigDev);
				cAxLim			= origBarHandle.Limits;
				
				
				axPredDev	= subplot(1,3,2,'Parent', panhandle,...
					'tag','axPredDev', 'NextPlot', 'add');
				contourDomainPlot(fem,domainID,yttU,1,axPredDev);
				perdBarHandle	=colorbar(axPredDev);
				perdBarHandle.Limits	=cAxLim;
				drawnow;
				
				axPredVar	= subplot(1,3,3,'Parent', panhandle,...
					'tag','axPredVar', 'NextPlot', 'add');
				contourDomainPlot(fem,domainID,VarYttU,1,axPredVar);
				
				%gif('partialMesSim.gif','frame',figComb,'DelayTime',1);
				drawnow;
			else
% 				figComb=plotStateAndVar(att,diag(ptt),rmseT(nSnaps),figComb);
				figComb=plotRmseErrorBarAndWrong(rmseT(1:nSnaps),sqrt(errorTVarDev(1:nSnaps)),perAllCorrect(1:nSnaps),figComb);

				cla(axPredDev);
				contourDomainPlot(fem,domainID,yttU,1,axPredDev);
				caxis(cAxLim)
				perdBarHandle	=colorbar(axPredDev);
				perdBarHandle.Limits	=cAxLim;
				drawnow;
				
				cla(axPredVar);
				contourDomainPlot(fem,domainID,VarYttU,1,axPredVar);
% 				hold all
				notCorrect	=find(~allCorrect);
				scatter3(nodeCoord(idMnp(notCorrect),1),nodeCoord(idMnp(notCorrect),2),nodeCoord(idMnp(notCorrect),3),'r','*')
				drawnow;
				wrongClassifData(:,1)	= devTest(partNum,idMnp(notCorrect));
				wrongClassifData(:,2)	= yttU(idMnp(notCorrect));
				wrongClassifData(:,3)	= pYttSave((notCorrect));
	 			wrongClassifData(:,4)	= sqrt(VarYttU(idMnp(notCorrect)));
				%gif
				
			end
			
			clear wrongClassifData;
			
		end		
	
		nSnaps=nSnaps+1;
		clear regionData;
	end	
	
	%% RMSE for (t+1) prediction
	if partNum < Nparts
		diffTp1				=bsxfun(@minus,devKey(partNum+1,:),ytp1);
	end
	madTp1				=mean(abs(diffTp1));
	mdTp1				=mean((diffTp1));
	rmseTp1(partNum)	=sqrt(sum((diff).^2)/size(devKey,2));  % rmse without noisy input
	
	seqPred(partNum).Snap(maxSnap).rmseTp1		=rmseTp1(partNum);
	seqPred(partNum).Snap(maxSnap).ErrorTp1Key	=diffTp1;
	seqPred(partNum).Snap(maxSnap).Ytp1			=ytp1U;
	seqPred(partNum).Snap(maxSnap).KeyPtID		=idMnp;
	seqPred(partNum).Snap(maxSnap).MadTp1		=madTp1;
	seqPred(partNum).Snap(maxSnap).MdTp1		=mdTp1;
	
	
	% predicting the state variables and their variances for (t+1)
	[attm1,pttm1]=getKalmanStatettm1(ptt,att,H,covErr);
	
end
%}

%%
% to plot the error contourplot for vnv case
% 
% for i=6 % choose pattern number
% 	
% 	cmin        =min(devTest(i,:));
% 	cmax        =max(devTest(i,:));
% 	cDelta      =0.05;     % step size for contourmap
% 	
% 	for kk=1:3% choose maxSnap
% 		
% 		if kk==maxSnap
% 			devST=devTest(i,:);
% 		else
% 			devST=seqPred(i).Snap(kk).Dev;
% 		end	
% 		
% 		allMesReg=seqPred(i).Snap(kk).Region;		
% 		madST(kk)=seqPred(i).Snap(kk).MadT;
% 		diffST	=devTest(i,:)-devST;		
% 		
% 		fig1        =figure;
% 		plotAxis1   =subplot(1,2,1);%gca;
% 		contourDomainPlot(fem,domainID,devTest(i,:),0,plotAxis1)
% 		view(0,0)                           % change figure viewpoint
% 		colormap(plotAxis1,'parula')
% % 		contourcmap('parula',cDelta);
% 		caxis([cmin cmax])
% 		
% % 		fig2        =figure;
% 		plotAxis1   =subplot(1,2,2);%gca;
% 		contourDomainPlot(fem,domainID,devST,0,plotAxis1)
% 		caxis([cmin cmax])
% 		view(0,0)   % change figure viewpoint
% 		colormap(plotAxis1,'parula')		
% % 		contourcmap('parula',cDelta);
% 		caxis([cmin cmax])
% 		
% 
% 		fig3        =figure;
% 		plotAxis1   =gca;
% 		contourDomainPlot(fem,domainID,diffST,0,plotAxis1)
% 		view(0,0)  
% 		contourcmap('jet',cDelta);
% 		cmin        =-4;%min(diffST);
% 		cmax        =4;%max(diffST);
% 		caxis([cmin cmax])
% 		colorbar
% 		
% 		
% 	end
% end