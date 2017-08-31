%% Kalman with krigging the error term, including functions for interpolating eigen values for inner
% AND CoMPARING THE effect on number of basis vectors with plotting the
% graphs in a separate section, also has provisions to check if predictions without kriging
% are better

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

sigmaMes    =1E-4;
nu.Type     ='diag';
nu.StdDev   =1e-4;
devAll      =devCenterd;
nBasis		=5:5:70;%5:5:50;    %     % set containing number of basis vectors


for k=1:length(nBasis)
	
	
	
	%% kalman recursion
	[keyEigVec,R,V,H,covErr]=getSystemMatricesSampled(nodeCoord,devAll,...
		sigmaMes,nu,iMnp,nBasis(k));
	devKey=devAll(:,iMnp);
	
	%% eigen interpolation
	fileNameCoarse= 'innerSelNodes - Copy.inp';
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
		[kt,att,ptt]=getKalmanStatett(pttm1,attm1,keyEigVec,R,V,zt(i,:)');
		% prediction
		[ytt(i,:),VarYtt(i,:),ytp1(i,:),VarYtp1(i,:)]=getPredictions(keyEigVec,att,ptt,H,V,R,covErr);
		
	end
	
	%% Adaptive and partial measurement
	tic
	seqPred(size(zt,1)).Snap(1).Dev=[];
	rmseTp1=zeros(size(zt,1),1);
	
	
	for i=numberCompleteMeasure+1:size(zt,1) % i for each part
		
		nSnaps=1;
		allMesReg=[];
		nRegions=1:size(nodeIDoutRectCoarse,2);
		rmseT=zeros(maxSnap,1);
		
		
		while ~isempty(nRegions)&&(nSnaps<=maxSnap)
			selIndex=zeros(length(nRegions),1);
			
			% prediction based on t-1 to predict first region to measure
			for j=1:length(nRegions)
				if nSnaps==1
					% when no region is measured use previous prediction (ytp1)
					devPart=ytp1(i-1,:);
					varPart=VarYtp1(i-1,:);
				else
					% when some regions are measured use updated prediction (ytt)
					devPart=ytt(i,:);
					varPart=VarYtt(i,:);
				end
				[devRegion,varRegion]=getRegionDevVarInner(nodeIDoutRectCoarse,nRegions(j),devPart,varPart,iMnp);
				% calculating the probability of being out of tolerance
				[selIndex(j)]=getRegionScore(tol, devRegion, varRegion);
			end
			
			[~,indexMax]=max(selIndex);
			selRegion=nRegions(indexMax);
			allMesReg=[allMesReg;selRegion]; % all the measured regions
			
			[tempV,tempR,tempEigVec,tempZ,regionIndex]=getRegionSysMatInner(nodeIDoutRectCoarse,...
				allMesReg,V,R,interpEigVec,devAll(i,:),iMnp);
			
			% updating the state variables after partial measurement
			[tempKt,att,ptt]=getKalmanStatett(pttm1,attm1,tempEigVec,tempR,tempV,tempZ');
			
			% predicting deviations for the whole part using updated state variables
			[ytt(i,:),VarYtt(i,:),ytp1(i,:),VarYtp1(i,:),yttU(i,:),ytp1U(i,:),VarYttU(i,:)]=getPredictionsNewInner(interpEigVec,...
				keyEigVec,tempEigVec,att,ptt,H,V,R,covErr,hyp,nodeCoord,...
				nodeCoord(iMnp,:),nodeCoord(regionIndex,:),tempZ);
% 			contourDomainPlot(fem,1,mpredictedU,1);
% 			% predicting deviations for the whole part using updated state
% 			% variables without krigging
% 			[ytt(i,:),VarYtt(i,:),ytp1(i,:),VarYtp1(i,:),yttU(i,:),ytp1U(i,:),VarYttU(i,:)]=getPredictionsNewInnerWoKrig(interpEigVec,...
% 				keyEigVec,tempEigVec,att,ptt,H,V,R,covErr,hyp,nodeCoord,...
% 				nodeCoord(iMnp,:),nodeCoord(regionIndex,:),tempZ);
			
			% removing the measured region from the list
			nRegions(indexMax)=[];
% 			RMSE calculation for the step
			diff=bsxfun(@minus,devCenterdKeypoint(i,:),ytt(i,:));
			rmseT(nSnaps)=sqrt(sum(diff.^2)/size(devCenterdKeypoint,2));  % rmse for key points
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
		diff=bsxfun(@minus,devCenterdKeypoint(i,:),ytp1(i-1,:));
		rmseTp1(i)=sqrt(sum((diff).^2)/size(devCenterdKeypoint,2));  % rmse without noisy input
		seqPred(i).Snap(maxSnap).rmseTp1=rmseTp1(i);
		seqPred(i).Snap(maxSnap).Ytp0=ytp1(i-1,:);
		
		% predicting the state variables and their variances for (t+1)
		[attm1,pttm1]=getKalmanStatettm1(ptt,att,H,covErr);
		
		
	end
	

	
	toc
	%% saving data
	fileString=sprintf('allDataInner_eigWithkrigYtp1%i.mat',nBasis(k));
	save(fileString,'seqPred');
% 	clear seqPred;
	
	%%  plots
	
% 		subPlotRos  =3;
% 		subPlotCol  =3;
% 		cmin        =-3;
% 		cmax        =3;
% 		nInstances  =1; % number of parts to plot
% 		cDelta=0.25;     % step size for contourmap
% 	
% 	
% 		for i=numberCompleteMeasure+5:numberCompleteMeasure+5+nInstances
% 	
% 			fig1        =figure;
% 			%     fig2        =figure;
% 			%     fig3        =figure;
% 			rmse=seqPred(i).Snap(maxSnap).RmseT;
% 			rmse=[rmseTp1(i);rmse];
% 	
% 	
% 			for kk=1:maxSnap
% 	
% 				devKey=seqPred(i).Snap(kk).Dev;
% 				var=seqPred(i).Snap(kk).Var;
% 				allMesReg=seqPred(i).Snap(kk).Region;
% 	
% 	
% 				% plot for prediction after each measurement
% 				set(0, 'currentfigure', fig1)
% 				plotAxis1=subplot(subPlotRos,subPlotCol,kk);
% 				contourDomainPlot(fem,idPart,devKey,0,plotAxis1)
% 				contourcmap('parula',cDelta);
% 				hold on
% 				allPoints=getCombinedregions(nodeIDoutRectCoarse,allMesReg);
% 				points2plot=nodeCoord(allPoints,:);
% 				scatter3(points2plot(:,1),points2plot(:,2),points2plot(:,3),'+k');
% 				view([0 0]);
% 				hold off
% 	
% 				%         % variance
% 				%         set(0, 'currentfigure', fig2)
% 				%         plotAxis2=subplot(subPlotRos,subPlotCol,kk);
% 				%         contourDomainPlot(fem,idPart,var,0,plotAxis2)
% 				%         %         contourcmap('parula',cDelta);
% 	
% 				%         %plot without measured region
% 				%         set(0, 'currentfigure', fig3)
% 				%         plotAxis1=subplot(subPlotRos,subPlotCol,kk);
% 				%         contourDomainPlot(fem,idPart,devKey,0,plotAxis1)
% 				%         contourcmap('parula',cDelta);
% 	
% 			end
% 	
% 	
% 			%     plotting original Deviation
% 			contourDomainPlot(fem,idPart,devAll(i,:),1)
% 			stringTitle=sprintf('Original variation pattern');
% 			title(stringTitle,'fontweight','bold','fontsize',13);
% 			hCont=contourcmap('parula',cDelta);
% 			%     hCont.FontSize=10;
% 	
% 			%plotting prediction for next time step
% 			contourDomainPlot(fem,idPart,ytp1U(i,:),1)
% 			stringTitle=sprintf('Prediction for next time step');
% 			title(stringTitle,'fontweight','bold','fontsize',13);
% 			hCont=contourcmap('parula',cDelta);
% 			%     hCont.FontSize=10;
% 	
% 	
% 			%plotting RMSE
% 			h=figure;
% 			myLine=plot(rmse);
% 			ax=h.Children;
% 			stringTitle=sprintf('RMSE for predicted deviations after each measurements');
% 			title(stringTitle,'fontweight','bold','fontsize',13,'FontName', 'AvantGarde');
% 			xlabel('Number of measurements','fontweight','bold');
% 			ylabel('RMSE in mm','fontweight','bold');
% 			set(ax, ...
% 				'Box'         , 'on'     , ...
% 				'xTick'       , 1:maxSnap+1 , ...
% 				'xTicklabel'  , 0:20, ...
% 				'yTick'       , 0:0.1:100 , ...
% 				'TickDir'     , 'in'     , ...
% 				'Ylim'        , [0 inf]     , ...
% 				'Xlim'        , [1 inf]     , ...
% 				'TickLength'  , [.01 .01] , ...
% 				'XMinorTick'  , 'off'     , ...
% 				'YMinorTick'  , 'off'     , ...
% 				'XColor'      , [.3 .3 .3], ...
% 				'YColor'      , [.3 .3 .3], ...
% 				'YGrid'       , 'on'      , ...
% 				'XGrid'       , 'on'      , ...
% 				'LineWidth'   , 1         ,...
% 				'FontName'    , 'Times New Roman',...
% 				'fontsize'    , 12         );
% 			myLine.LineWidth=1.5;
% 			myLine.Color='b' ;
% 	
% 	
% 			fig1=makeSubplotCompact(fig1,subPlotRos, subPlotCol);
% 			fig1 = tightfig(fig1);
% 			set(0, 'currentfigure', fig1)
% 			stringTitle=sprintf('Sequence of adaptive measurements based on information gain');
% 			figtitle(stringTitle);
% 	
% % 			fig2=makeSubplotCompact(fig2,subPlotRos, subPlotCol);
% % 			fig2 = tightfig(fig2);
% % 			set(0, 'currentfigure', fig2)
% % 			stringTitle=sprintf('Variance of predicted deviations after each measurement');
% % 			figtitle(stringTitle);
% % 
% % 			fig3=makeSubplotCompact(fig3,subPlotRos, subPlotCol);
% % 			fig3 = tightfig(fig3);
% % 			set(0, 'currentfigure', fig3)
% % 			stringTitle=sprintf('Predicted complete part deviations after each measurement');
% % 			figtitle(stringTitle);
% % 
% % 			% for single colorbar
% % 			set(0, 'currentfigure', fig2)
% % 
% % 			h=colorbar;
% % 			h.FontSize=12;
% % 			set(get(h,'Ylabel'),'string','mm^{2}','fontweight','bold','fontsize',14)
% % 			set(h, 'Position', [.8314 .11 .03 .8]);
% % 
% % 			set(0, 'currentfigure', fig3)
% % 			h=colorbar;
% % 			h.FontSize=12;
% % 			set(get(h,'Ylabel'),'string','mm','fontweight','bold','fontsize',14)
% % 			set(h, 'Position', [.8314 .11 .03 .8])
% 			
% 	
% 		end
	
%  Plotting eigen basis
		
% 		for i=1:size(keyEigVec,2)
% 		    contourDomainPlot(fem,idPart,keyEigVec(:,i),1);
% 		    print(sprintf('Basis%i',i),'-dpng')
% 		end
% 	
% 		% %Plotting regions and key points
% 	
% 		meshplotAxisDefined(fem,1);
% 		hold all;
% 		figDense=gcf;
% 	
% 		meshplotAxisDefined(fem,1);
% 		hold all;
% 		figCoarse=gcf;
% 		for i=1:size(nodeIDoutRectCoarse,2)
% 	
% 			set(0, 'currentfigure', figDense)
% 			regionIndex=find(nodeIDoutRectCoarse(:,i));%% finds selected nodes
% 			scatter3(nodeCoord(regionIndex,1),nodeCoord(regionIndex,2),nodeCoord(regionIndex,3),...
% 				'o','filled','linewidth',2);
% 	
% 			set(0, 'currentfigure', figCoarse)
% 			regionIndex=find(nodeIDoutRectDense(:,i));%% finds selected nodes
% 			scatter3(nodeCoord(regionIndex,1),nodeCoord(regionIndex,2),nodeCoord(regionIndex,3),...
% 				'o','filled','linewidth',2);
% 	
% 		end


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
for k=1:length(nBasis)
	
	fileString=sprintf('allDataInner_eigWithkrigYtp1%i.mat',nBasis(k));
	seqPred=load (fileString);
	seqPred=seqPred.seqPred;
	
	for i=numberCompleteMeasure+1:length(seqPred)
		rmsePlot(:,i)=seqPred(i).Snap(maxSnap).RmseT;		% accessing each of the 50 replication
	end
	
	rmseAvg(:,k)=sum(rmsePlot,2)./40; % averaging all replications
	
end

% plotting eig cmparison
rsmeFig1=figure;
lineRmseEig=plot(rmseAvg);

legRmseEig=cell(length(lineRmseEig),1);


for k=1:length(lineRmseEig)
	lineRmseEig(k).LineWidth=1.5;
	lineRmseEig(k).Color=cmapRmsePlot(k,:);
	lineRmseEig(k).Marker=markers{k};
	lineRmseEig(k).MarkerFaceColor=cmapRmsePlot(k,:);
% 	lineRmseEig(k).MarkerEdgeColor='none';
	legRmseEig{k}=sprintf('%i',nBasis(k));
	
end

stringTitle='Average RMSE comparison after each measurement for different number of basis';
stringXlabel='Number of measurements';
stringYlabel='RMSE in mm';
rsmeFig1=changeAxesLooks(rsmeFig1,'',stringXlabel,stringYlabel);
[~] = legendflex(lineRmseEig, legRmseEig, ...% 'ref', hl(2).leg, ...
	'anchor'    ,{'ne','ne'}, ...
	'buffer'    ,[-10 -5]      , ...
	'nrow'      ,4          ,...
	'fontsize'  ,12         , ...
	'box'       ,'off'       ,...
	'title'     ,'Number of basis shapes'   );
% rsmeFig1('PaperPositionMode','auto');
print(stringTitle,'-dpdf');

% plotting absolute eig comparison
fig=figure;
rmseAvgeEig=sum(rmseAvg)./20;  % summming errors from all partial measurements for a given number of eigen vectors
lineRmseEigAvg=plot(rmseAvgeEig);
fig=changeAxesLooks(fig,'Average RMSE for different number of eign vectors ',...
    'Number of eigen basis','RMSE in mm');
lineRmseEigAvg.LineWidth=1.5;
plotAxis=fig.Children;
plotAxis.XTickLabel=num2cell(nBasis);
% fig('PaperPositionMode','auto');
print('eigCompInnerAvg','-dpdf');

% % plotting absolute first five measurements
% fig=figure;
% rmseAvg(6:end,:)=[];
% rmseAvgeEig=sum(rmseAvg)./20;  % summming errors from all partial measurements for a given number of eigen vectors
% lineRmseEigAvg=plot(rmseAvgeEig);
% fig=changeAxesLooks(fig,'Average RMSE for different number of eign vectors ',...
%     'Number of eigen basis','RMSE in mm');
% lineRmseEigAvg.LineWidth=1.5;
% plotAxis=fig.Children;
% plotAxis.XTickLabel=num2cell(nBasis);

%% plotting zeroth prediction and priginal pattern

patternNo=34;
cDelta= 0.25;
% cmin        =1.1*min(devPatterns(patternNo,:));
% cmax        =1.1*max(devPatterns(patternNo,:));
devST0=seqPred(patternNo).Snap(maxSnap).Ytp0;
contourDomainPlot(fem,idPart,devST0,1)
contourcmap('parula',cDelta);
% caxis([cmin cmax])
view(0,0)
print('zerothPred','-dpdf');

devST=devPatterns(patternNo,:);
contourDomainPlot(fem,idPart,devST,1)
view(0,0)
contourcmap('parula',cDelta);
% fig('PaperPositionMode','auto');
print('originalDev','-dpdf');

%% debug code
% for test=11:50
% 	diff(test,:)=bsxfun(@minus,devCenterdKeypoint(test,:),ytt(test,:));
% 	rmseT(nSnaps)=sqrt(sum(diff(test,:).^2)/size(devCenterdKeypoint,2));  % rmse for key points
% end
