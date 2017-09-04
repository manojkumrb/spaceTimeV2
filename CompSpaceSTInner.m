% code for comparision of space time KF and krigging


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
partRegions				=load('inner9Regions.mat');%inner_regions.mat');% contains variable nodeIDoutRect
nodeIDoutRect		=partRegions.nodeIDoutRect;


devPatterns				=load('simAutoCorDevInnerBatchesCombined.mat');%  %hingeDevArSimLoc
devPatterns				=devPatterns.simData;

seqPred					=load('allDataInner_eigWithkrigYtp140.mat');%predEigBatch2NumBasis30WithKrigBatchComb.mat');
seqPred					=seqPred.seqPred;

%% loading from predefined selection
selNodes		=load('doorInnerSelNodes.mat');
selNodes		=selNodes.selNodes;
Mnp				=selNodes.Coord;
iMnp			=selNodes.ID;

% devCenterd			=devPatterns;%devCenterd1+randn(size(devPatterns))*10e-3;%devPatterns;x=z*sigma+mu;; % set equal to input in a coarse mesh
% devCenterdKeypoint  =devCenterd(:,iMnp);
%% setting test and train data
nTrainBatch			=[1,2,3,4,5,6,7,8];    % set of batches to run with each value represnting number of consecutive batches to add in train data

%% learning parameters
%%%%%%%%%%%%%%%%%%%%%%%%%
sourceBatch		=6;			% coming from the batch used to test predictions in batch size comparison
sourcePattern   =18;%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
x               = nodeCoord; % create a data set
devSmooth		=devPatterns(sourceBatch).FemDevDomain;

%% mean and covariance setting
m0              ={'meanZero'};  hyp.mean = [];                      % no hyperparameters are needed
sn              =0.001;                                             % noise standard dev
likfunc         ={'likGauss'};   hyp.lik = log(sn);                 % gaussian likelihood
sf              =1;                                                 % latent function standard dev
lScale          =[60;50;100]*rand(1);                               %rand(D,1);% characteristic length scale
covMatern       ={@covMaternard,5};  hypMat=log([lScale;sf]);      % Matern class d=5
% covMatern       ={@covSEard};  hypMat=log([lScale;sf]);      % modified matern to sq exp
inf				= @infGaussLik;

hyp.cov         =hypMat;

xLimited        =Mnp;% sampleNodesDomain(fem,domainID,nodeLimit);

cov = {'apxSparse',covMatern,xLimited};           % inducing points
infv  = @(varargin) inf(varargin{:},struct('s',01.0));           % VFE, opt.s = 0; % FITC, opt.s = 1;% SPEP, 0<opt.s<1


%% log marginal likelihood optimization
% fitc is used here to avoid the trouble of finding key points again
hyp				= load('hypInnerFitc.mat');%minimize(hyp, @gp, -50,  infv, m0, cov, likfunc,  x, devSmooth(sourcePattern,:)'); % for covFITC
hyp             = hyp.hyp;
%% running spatial gpr prediction
numberCompleteMeasure   =10;
maxSnap                 =length(seqPred(end).Snap);
nPattern                =size(devSmooth,1);
gprPred(nPattern).Snap(1).Dev=[];


% % loading gpr predictions
% load('gprPredSourceBatc6InnerMatenCov.mat');
% load('gprPredSourceBatc6InnerSEcov.mat');


for i=numberCompleteMeasure+1:nPattern % numberCompleteMeasure+5 % i for each part
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    inputPattern=i;%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    rmse=zeros(maxSnap,1);
    
    for j=1:maxSnap
        
        %% Selecting region with measurement
        allMesReg       =seqPred(inputPattern).Snap(j).Region;  % gets regions measured upto a given number of snap
        allPoints       =zeros(size(nodeIDoutRect,1),1);
        for k =1:length(allMesReg)
            allPoints   =allPoints|nodeIDoutRect(:,allMesReg(k)); % aggregates all key points from measured regions
        end
        %%%Plot to verify region selection
        %         meshplotAxisDefined(fem,domainID )
        %         hold all
        %         plotID=find(allPoints);
        %         plot3(nodeCoord(plotID,1),nodeCoord(plotID,2),nodeCoord(plotID,3),'*')
        %         hold off
        
        xt=nodeCoord(allPoints,:);
        yt=devSmooth(inputPattern,allPoints);
        % [xt,sample]=sampleNodes(xt,0.9);
        % yt=yt(sample);
        
        z=nodeCoord;
        hypPred.cov=hyp.cov;
        hypPred.lik=hyp.lik;
		[mpredicted, s2] = gp(hypPred,infv,m0,cov,likfunc,xt, yt', z);  % with cov funtion modified for fitc

%         [mpredicted, s2] = gp(hypPred, @infExact, m0, csu, likfunc, xt, yt', z);
        mpredicted=mpredicted';
        
        difference=devSmooth(inputPattern,:)-mpredicted;
        rmse(j)=sqrt(sum(difference.^2)/length(difference));
        gprPred(inputPattern).Snap(j).Rmse=rmse;
        gprPred(inputPattern).Snap(j).Dev=mpredicted;
        
        
    end
    
end


%% plotting
subPlotRos  =1;
subPlotCol  =3;
nInstances  =30; % number of parts to plot


for i=numberCompleteMeasure+24%numberCompleteMeasure+1:numberCompleteMeasure+1+nInstances %
    cmin        =1.1*min(devSmooth(i,:));
    cmax        =1.1*max(devSmooth(i,:));
    cDelta      =0.25;     % step size for contourmap
    
    rmseST=seqPred(i).Snap(maxSnap).RmseT;
    rmseSTp1=seqPred(i).Snap(maxSnap).rmseTp1;    % getting rmse with no measurement
    rmseST=[rmseSTp1;rmseST];
    rmseS=gprPred(i).Snap(maxSnap).Rmse;
    rmseS=[NaN;rmseS];
    
%% % to plot the contourplot comparison between space and spacetime    
%     for kk=1:maxSnap
%         
%         devST=seqPred(i).Snap(kk).Dev;
%         devS=gprPred(i).Snap(kk).Dev;
%         
%         var=seqPred(i).Snap(kk).Var;
%         allMesReg=seqPred(i).Snap(kk).Region;
%         
%         fig1        =figure;
% %         % plot for prediction after each measurement
%         set(0, 'currentfigure', fig1)
%         plotAxis1   =subplot(subPlotRos,subPlotCol,1);
%         contourDomainPlot(fem,idPart,devSmooth(i,:),0,plotAxis1)
%         caxis([cmin cmax])
%         view(0,0)                           % change figure viewpoint
%         
%         plotAxis1   =subplot(subPlotRos,subPlotCol,2);
%         contourDomainPlot(fem,idPart,devS,0,plotAxis1)
%         caxis([cmin cmax])
%         view(0,0)                           % change figure viewpoint
%         
%         plotAxis1   =subplot(subPlotRos,subPlotCol,3);
%         contourDomainPlot(fem,idPart,devST,0,plotAxis1)
%         caxis([cmin cmax])
%         view(0,0)                           % change figure viewpoint
%         
%         contourcmap('parula',cDelta);
%         fig1=makeSubplotCompact(fig1,subPlotRos, subPlotCol);
%         
% % %       saving figure
%         figname=sprintf('innerSpaceStCompSnap%i.png',kk);
% %         print('-r400',fig1,figname,'-dpng')
%     end
    
%% Plot RMSE for individual pattern
    h=figure;
    myLine2=plot(rmseST);
    hold all;
    myLine=plot(1:length(rmseST),rmseS);
    
    
    
    stringTitle=sprintf('RMSE for predicted deviations after each measurements');
    stringXlabel='Number of measurements';
    stringYlabel='RMSE in mm';
    h=changeAxesLooks(h,stringTitle,stringXlabel,stringYlabel);    
    myLine.LineWidth=1.5;
    myLine2.LineWidth=1.5;
    % legend
    legendflex([myLine2,myLine], {'Spatio-temporal Kalman filter','Spatial Kriging'}, ...
        'anchor', {'ne','ne'}, ...
        'buffer', [0 0], ...
        'nrow',   2,...
        'fontsize',10, ...
        'xscale', 0.66, ...
        'box','on');
	
	fileName=sprintf('spaceAndSTCompInnerPttern%i',i);
	print('-r400',fileName,'-dpng');
end
%% Plot to select a colourful individual pattern
% 
% for i=numberCompleteMeasure+1:nPattern
% 	contourDomainPlot(fem,idPart,devSmooth(i,:),1);
% end
% 
% 


%% Rmse diff AFTER EACH MEASUREMET AVERAGED FOR ALL PATTERNS

rmseST=zeros(1,maxSnap);
rmseS=zeros(1,maxSnap);

for kk=1:maxSnap
    diff=0;

    for i=numberCompleteMeasure+1:nPattern
        
        rmseST(kk)=rmseST(kk)+seqPred(i).Snap(kk).RmseT(kk);
        rmseS(kk)=rmseS(kk)+gprPred(i).Snap(kk).Rmse(kk);
        
    end
           
end

totalInstances= nPattern-(numberCompleteMeasure+1);

delta=((rmseS-rmseST)./rmseS)*100;
fig=figure;
hh=plot(delta);
fig=changeAxesLooks(fig,'','Number of measurements','Percentage improvement');
print('spaceAndSTCompInnerImprovePer','-dpdf');

fig=figure;
hh=plot(1:length(rmseS),rmseS./totalInstances,1:length(rmseST),rmseST./totalInstances);
fig=changeAxesLooks(fig,'','Number of measurements','RMSE in mm');

set(hh,'linewidth',1.5);

legendflex(hh, {'Spatial prediction','Spatio-temporal prediction'}, ...
        'anchor', {'ne','ne'}, ...
        'buffer', [0 0], ...
        'nrow',   2,...
        'fontsize',12, ...
        'xscale', 0.66, ...
        'box','on');
print('spaceAndSTCompInner','-dpdf');

