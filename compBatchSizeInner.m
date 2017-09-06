% modifications for checking the effect of learning from more batches
close all;
clear ;
% dbstop if error;

run 'C:\Users\babu_m\Documents\GitHub\source\initFiles.m';
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

devPatterns     =load('simAutoCorDevInnerBatchesCombined.mat');%  %hingeDevArSimLoc
devPatterns     =devPatterns.simData;

%% loading from predefined selection
selNodes=load('doorInnerSelNodes.mat');
selNodes=selNodes.selNodes;
Mnp=selNodes.Coord;
iMnp=selNodes.ID;
devTest				= devPatterns(8).FemDevDomain;   % fixed to 6 for now
devCenterdKeypoint  = devTest(:,iMnp);
%% setting test and train data
nTrainBatch			=[1,2,3,4,5,6,7];    % set of batches to run with each value represnting number of consecutive batches to add in train data
% stopped at 8th batch, as it is used as test case

for j=1:length(nTrainBatch)       % batches to run
    
    
    
    devTrain			=[];
    
    for l=1:nTrainBatch(j)        % adding train data
        devTrain		=[devTrain;devPatterns(l).FemDevDomain];
    end
    
   
    
    %% finding key points in each region
    keyPointIndex		  = zeros(size(nodeIDoutRectDense,1),1);
    keyPointIndex(iMnp)	  = 1;
    nodeIDoutRectCoarse= bsxfun(@times,keyPointIndex,nodeIDoutRectDense);
    
    %% load optimised error process covariance parameters
    load('hypInnerFitc.mat');%optimisedHypParmHinge.mat'); %contains struct hyp in gp function format(i.e. hyp.cov=log(parameters)) for matern covariance
    
    sigmaMes    =1E-4;
    nu.Type     ='diag';
    nu.StdDev   =1e-4;
    nBasis		=35;%[2,3,4,5,6,7,8,9,10:5:60];         % set containing number of basis vectors to be tested
    
    
    for k=1:length(nBasis)
        
        sprintf('running %i batch %i basis',j,nBasis(k))
        
        %% kalman recursion
        [keyEigVec,R,V,H,covErr]=getSystemMatricesSampled(nodeCoord,devTrain,...
            sigmaMes,nu,iMnp,nBasis(k));
        
        
        
        %% eigen interpolation
        fileNameCoarse= 'innerSelNodes - Copy.inp';
        interpEigVec=getEigenInterp(nodeCoord,fileNameCoarse,keyEigVec, domainID);
        % interpEigVec=load('interpEigVecInner3_400.mat');
        % interpEigVec=interpEigVec.interpEigVec;
        %%
        devKey		=devTest(:,iMnp);
        zt          =devKey; % measurements
        tol         =4;    % is a tricky value affects the adaptivity a lot..currently works well when tol is set to 2; set to for after sensitivity for inner 9 regions
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
            [kt,att,ptt]=getKalmanStatett(pttm1,attm1,keyEigVec,R,V,zt(i,:)');
            % prediction
            [ytt(i,:),VarYtt(i,:),ytp1(i,:),VarYtp1(i,:)]=getPredictions(keyEigVec,att,ptt,H,V,R,covErr);
            
        end
        
        %% Adaptive and partial measurement
        tic
        
        seqPred(size(zt,1)).Snap(1).Dev=[];
        rmseTp0=zeros(size(zt,1),1);
        
        
        for i=numberCompleteMeasure+1:size(zt,1) % i for each part
            
            nSnaps=1;
            allMesReg=[];
            nRegions=1:size(nodeIDoutRectCoarse,2);
            rmseT=zeros(maxSnap,1);
            
            
            while ~isempty(nRegions)&&(nSnaps<=maxSnap)
                selIndex=zeros(length(nRegions),1);
                
                % prediction based on t-1 to predict first region to measure
                for jj=1:length(nRegions)
                    if nSnaps==1
                        % when no region is measured use previous prediction (ytp1)
                        devPart=ytp1(i-1,:);
                        varPart=VarYtp1(i-1,:);
                    else
                        % when some regions are measured use updated prediction (ytt)
                        devPart=ytt(i,:);
                        varPart=VarYtt(i,:);
                    end
                    [devRegion,varRegion]=getRegionDevVarInner(nodeIDoutRectCoarse,nRegions(jj),devPart,varPart,iMnp);
                    % calculating the probability of being out of tolerance
                    [selIndex(jj)]=getRegionScore(tol, devRegion, varRegion);
                end
                
                [~,indexMax]=max(selIndex);
                selRegion=nRegions(indexMax);
                allMesReg=[allMesReg;selRegion]; % all the measured regions
                
                [tempV,tempR,tempEigVec,tempZ,regionIndex]=getRegionSysMatInner(nodeIDoutRectCoarse,...
                    allMesReg,V,R,interpEigVec,devTest(i,:),iMnp);
                
                % updating the state variables after partial measurement
                [tempKt,att,ptt]=getKalmanStatett(pttm1,attm1,tempEigVec,tempR,tempV,tempZ');
                
                % % predicting deviations for the whole part using updated
                %state variables with krigging error included
                [ytt(i,:),VarYtt(i,:),ytp1(i,:),VarYtp1(i,:),yttU(i,:),ytp1U(i,:),VarYttU(i,:)]=getPredictionsNewInner(interpEigVec,...
                    keyEigVec,tempEigVec,att,ptt,H,V,R,covErr,hyp,nodeCoord,...
                    nodeCoord(iMnp,:),nodeCoord(regionIndex,:),tempZ);

                % removing the measured region from the list
                nRegions(indexMax)=[];
                % RMSE calculation for the step
				diff=bsxfun(@minus,devCenterdKeypoint(i,:),ytt(i,:));
				rmseT(nSnaps)=sqrt(sum((diff).^2)/size(devCenterdKeypoint,2));  % rmse without noisy input

                
                % saving sequential predictions
                seqPred(i).Snap(nSnaps).Dev=yttU(i,:);
                seqPred(i).Snap(nSnaps).Var=VarYttU(i,:);
                seqPred(i).Snap(nSnaps).RmseT=rmseT;
                seqPred(i).Snap(nSnaps).Region=allMesReg;%selRegion;
                nSnaps=nSnaps+1;
            end
            
            % RMSE for (t+1) prediction
            diff=bsxfun(@minus,devCenterdKeypoint(i,:),ytp1(i-1,:));
%             diff=bsxfun(@minus,devTest(i,:),ytp1U(i-1,:));          % rmse with zero measurements on the ith pattern
            rmseTp0(i)=sqrt(sum((diff).^2)/size(devTest,2));        % rmse without noisy input
            seqPred(i).Snap(maxSnap).rmseTp0=rmseTp0(i);
            
            % predicting the state variables and their variances for (t+1)
            [attm1,pttm1]=getKalmanStatettm1(ptt,att,H,covErr);
            
            
        end
        
        toc
        % saving data
        fileString=sprintf('predEig%iBatch%iBasisWithKrig.mat',nTrainBatch(j) ,nBasis(k));
        save(fileString,'seqPred','-v7.3');
        
    end
    
end



%%  reading data from the batch analysis files and plotting

nBatches	=[1,2,3,4,5,6,7];
% plot parameters
markers = {'o','s','d','^','v','x','+','*','.','>','<','p','h','o','s','d','^','v','x'};
cmapRmsePlot = cbrewer('qual','Set1',length(nBatches));

rmseAvg=[];
rmsePlot=[];
nInstances=length(seqPred)-(numberCompleteMeasure+1);
for k=1:length(nBatches)
    
    fileString=sprintf('predEig%iBatch35BasisWithKrig.mat',nBatches(k));
    seqPred=load (fileString);
    seqPred=seqPred.seqPred;
    
    for i=numberCompleteMeasure+1:length(seqPred)
        rmsePlot(:,i)=seqPred(i).Snap(maxSnap).RmseT;		% accessing each of the 50 replication
    end
    
    rmseAvg(:,k)=sum(rmsePlot,2)./nInstances; % averaging all replications
    
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
    legRmseEig{k}=sprintf('%i',nBatches(k));
    
end

stringTitle=[{'RMSE comparison after each measurement for'},{'different number of batches learnt from'}];
stringXlabel='Number of measurements';
stringYlabel='RMSE in mm';
rsmeFig1=changeAxesLooks(rsmeFig1,stringTitle,stringXlabel,stringYlabel);
[~] = legendflex(lineRmseEig, legRmseEig, ...% 'ref', hl(2).leg, ...
    'anchor'    ,{'ne','ne'}, ...
    'buffer'    ,[-10 -5]      , ...
    'nrow'      ,2          ,...
    'fontsize'  ,12         , ...
    'box'       ,'off'       ,...
    'title'     ,'Number of batches learnt from'   );
print('-r400','batchSizeCompInnerAvg','-dpng');


% plotting absolute eig comparison
fig=figure;
rmseAvgeEig=sum(rmseAvg)./20;  % summming errors from all partial measurements for a given number of eigen vectors
lineRmseEigAvg=plot(rmseAvgeEig);
fig=changeAxesLooks(fig,'Average RMSE versus different number of batches learnt from ',...
    'Number of batches learnt from','RMSE in mm');
lineRmseEigAvg.LineWidth=1.5;
plotAxis=fig.Children;
plotAxis.XTickLabel=num2cell(nBatches);
print('-r400','batchSizeSensitInnerAvg','-dpng');
