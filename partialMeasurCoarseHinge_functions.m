%%Partial measurement data case implementaion for hinge, tidy function
%%implementaion
close all;
clear ;
% dbstop if error;

run 'C:\Users\babu_m\Documents\GitHub\source\initFiles.m';
eps         =10E-16;
fem         =femInit();
fileName{1} =strcat(inpFiles,'\Remesh_hinge.inp');%\top_hat.inp');%\halo.inp');%
% fem       =importMesh(fem,fileName{1});%'top_hat.inp'); %use only for top hat
fem         =importMultiMesh(fem, fileName);%has provisions for domain separation within code
fem         =femPreProcessing(fem);
domainID    =1;
idPart      =domainID;


coordi          =fem.xMesh.Node.Coordinate;
nodeIdDomain    =fem.Domain(idPart).Node;
nodeCoord       =fem.xMesh.Node.Coordinate(nodeIdDomain,:); % coordinates for the domain
% load pre-defined part regions
partRegions     =load('hingeRegionsCoarse.mat');% contains variable nodeIDoutRect
nodeIDoutRect   =partRegions.nodeIDoutRect;
% load optimised error process covariance parameters
load('optimisedHypParm.mat'); %contains struct hyp in gp function format(i.e. hyp.cov=log(parameters)) for matern covariance
% for coarse mesh
devPatterns     =load('hingeDevArSimLocCoarse.mat');%  %hingeDevArSimLoc
devPatterns     =devPatterns.femDevDomain;
% % % testing the effect of not centering
devCenterd1     =devPatterns;
% devCenterd1=bsxfun(@minus,devPatterns,mean(devPatterns)); % to calculate the rmse and for modelling deviation from mean
% devSelected= devPatterns(:,iMnp);
devCenterd      =devCenterd1;%devCenterd1+randn(size(devPatterns))*10e-3;%devPatterns;x=z*sigma+mu;; % set equal to input in a coarse mesh

sigmaMes    =1E-4;
nu.Type     ='diag';
nu.StdDev   =1e-4;
dev         =devCenterd;

%% kalman recursion
[eigVec,R,V,H,covErr]=getSystemMatrices(nodeCoord,dev,sigmaMes,nu);

zt          =devCenterd; % measurements
tol         =2;    % is a tricky value affects the adaptivity a lot..currently works well when tol is set to 2
countPart   =1;
maxSnap     =6;
numberCompleteMeasure= 10; % to stabalize the prediction

pc          =dev*eigVec;
% at time zero
a00         =pc(1,:)';
p00         =eye(size(eigVec,2));
ytp1        =zeros(size(devCenterd,1),size(nodeCoord,1));
ytt         =zeros(size(ytp1));
VarYtt      =zeros(size(ytp1));
VarYtp1     =zeros(size(ytp1));


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
    [kt,att,ptt]=getKalmanStatett(pttm1,attm1,eigVec,R,V,zt(i,:)');
    % prediction
    [ytt(i,:),VarYtt(i,:),ytp1(i,:),VarYtp1(i,:)]=getPredictions(eigVec,att,ptt,H,V,R,covErr);
    
end

%% Adaptive and partial measurement

seqPred(size(zt,1)).Snap(1).Dev=[];
rmseTp1=zeros(size(zt,1),1);


for i=numberCompleteMeasure+1:size(zt,1) % i for each part
    
    nSnaps=1;
    allMesReg=[];
    nRegions=1:size(nodeIDoutRect,2);
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
            [devRegion,varRegion]=getRegionDevVar(nodeIDoutRect,nRegions(j),devPart,varPart);
            % calculating the probability of being out of tolerance
            [selIndex(j)]=getRegionScore(tol, devRegion, varRegion);
        end
        
        [~,indexMax]=max(selIndex);
        selRegion=nRegions(indexMax);
        allMesReg=[allMesReg;selRegion]; % all the measured regions
        [tempV,tempR,tempEigVec,tempZ,regionIndex]=getRegionSysMat(nodeIDoutRect,allMesReg,V,R,eigVec,zt(i,:));
        % updating the state variables after partial measurement
        [tempKt,att,ptt]=getKalmanStatett(pttm1,attm1,tempEigVec,tempR,tempV,tempZ');
        % predicting deviations for the whole part using updated state variables
        [ytt(i,:),VarYtt(i,:),ytp1(i,:),VarYtp1(i,:)]=getPredictions(eigVec,att,ptt,H,V,R,covErr);
        
        % removing the measured region from the list
        nRegions(indexMax)=[];
        % RMSE calculation for the step
        diff=bsxfun(@minus,devCenterd1(i,:),ytt(i,:));
        rmseT(nSnaps)=sqrt(sum((diff).^2)/size(devCenterd1,2));  % rmse without noisy input
        
        % saving sequential predictions
        seqPred(i).Snap(nSnaps).Dev=ytt(i,:);
        seqPred(i).Snap(nSnaps).Var=VarYtt(i,:);
        seqPred(i).Snap(nSnaps).RmseT=rmseT;
        seqPred(i).Snap(nSnaps).Region=allMesReg;%selRegion;
        nSnaps=nSnaps+1;
    end
    
    % RMSE for (t+1) prediction
    diff=bsxfun(@minus,devCenterd1(i,:),ytp1(i-1,:));
    rmseTp1(i)=sqrt(sum((diff).^2)/size(devCenterd1,2));  % rmse without noisy input
    seqPred(i).Snap(maxSnap).rmseTp1=rmseTp1(i);
    
    % predicting the state variables and their variances for (t+1)
    [attm1,pttm1]=getKalmanStatettm1(ptt,att,H,covErr);
    
    
end


%% new plots

subPlotRos  =2;
subPlotCol  =3;
cmin        =-3;
cmax        =3;
nInstances  =3; % number of parts to plot
cDelta=0.25;     % step size for contourmap


for i=numberCompleteMeasure+1:numberCompleteMeasure+nInstances
    
    fig1        =figure;
    fig2        =figure;
    fig3        =figure;
    rmse=seqPred(i).Snap(maxSnap).RmseT;
    rmse=[rmseTp1(i);rmse];
    
    
    for kk=1:maxSnap
        
        dev=seqPred(i).Snap(kk).Dev;
        var=seqPred(i).Snap(kk).Var;
        allMesReg=seqPred(i).Snap(kk).Region;
        
        
        % plot for prediction after each measurement
        set(0, 'currentfigure', fig1)
        plotAxis1=subplot(subPlotRos,subPlotCol,kk);
        contourDomainPlot(fem,idPart,dev,0,plotAxis1)
        contourcmap('parula',cDelta);
        hold on
        allPoints=getCombinedregions(nodeIDoutRect,allMesReg);
        points2plot=nodeCoord(allPoints,:);
        scatter3(points2plot(:,1),points2plot(:,2),points2plot(:,3),'+k');
        hold off
        
        set(0, 'currentfigure', fig2)
        plotAxis2=subplot(subPlotRos,subPlotCol,kk);
        contourDomainPlot(fem,idPart,var,0,plotAxis2)
        %         contourcmap('parula',cDelta);
        %plot without measured region
        set(0, 'currentfigure', fig3)
        plotAxis1=subplot(subPlotRos,subPlotCol,kk);
        contourDomainPlot(fem,idPart,dev,0,plotAxis1)
        contourcmap('parula',cDelta);
        
    end
    
    
    %plotting original Deviation
    contourDomainPlot(fem,idPart,zt(i,:),1)
    stringTitle=sprintf('Original variation pattern');
    title(stringTitle,'fontweight','bold','fontsize',13);
    hCont=contourcmap('parula',cDelta);
%     hCont.FontSize=10;
    
    %plotting prediction for next time step
    contourDomainPlot(fem,idPart,ytp1(i,:),1)
    stringTitle=sprintf('Prediction for next time step');
    title(stringTitle,'fontweight','bold','fontsize',13);
    hCont=contourcmap('parula',cDelta);
%     hCont.FontSize=10;
    
    
    %plotting RMSE
    h=figure;
    myLine=plot(rmse);
    ax=h.Children;
    stringTitle=sprintf('RMSE for predicted deviations after each measurements');
%     title(stringTitle,'fontweight','bold','fontsize',13,'FontName', 'AvantGarde');
    xlabel('Number of measurements','fontweight','bold');
    ylabel('RMSE in mm','fontweight','bold');
    set(ax, ...
        'Box'         , 'on'     , ...
        'xTick'       , 1:maxSnap+1 , ...
        'xTicklabel'  , ['0';'1';'2';'3';'4'] , ...
        'yTick'       , 0:0.1:100 , ...
        'TickDir'     , 'in'     , ...
        'Ylim'        , [0 inf]     , ...
        'Xlim'        , [1 inf]     , ...
        'TickLength'  , [.01 .01] , ...
        'XMinorTick'  , 'off'     , ...
        'YMinorTick'  , 'off'     , ...
        'XColor'      , [.3 .3 .3], ...
        'YColor'      , [.3 .3 .3], ...
        'YGrid'       , 'on'      , ...
        'XGrid'       , 'on'      , ...
        'LineWidth'   , 1         ,...
        'FontName'    , 'AvantGarde',...
        'fontsize'    , 12         );
    myLine.LineWidth=1.5;
    myLine.Color='b' ;
    
    
    fig1=makeSubplotCompact(fig1,subPlotRos, subPlotCol);
    fig1 = tightfig(fig1);
    set(0, 'currentfigure', fig1)
    stringTitle=sprintf('Sequence of adaptive measurements based on information gain');
    figtitle(stringTitle);
    
    fig2=makeSubplotCompact(fig2,subPlotRos, subPlotCol);
    fig2 = tightfig(fig2);
    set(0, 'currentfigure', fig2)
    stringTitle=sprintf('Variance of predicted deviations after each measurement');
    figtitle(stringTitle);
    
    fig3=makeSubplotCompact(fig3,subPlotRos, subPlotCol);
    fig3 = tightfig(fig3);
    set(0, 'currentfigure', fig3)
    stringTitle=sprintf('Predicted complete part deviations after each measurement');
    figtitle(stringTitle);
    
    % for single colorbar
    set(0, 'currentfigure', fig2)
    
    h=colorbar;
    h.FontSize=12;
    set(get(h,'Ylabel'),'string','mm^{2}','fontweight','bold','fontsize',14)
    set(h, 'Position', [.8314 .11 .03 .8]);
    
    set(0, 'currentfigure', fig3)
    h=colorbar;
    h.FontSize=12;
    set(get(h,'Ylabel'),'string','mm','fontweight','bold','fontsize',14)
    set(h, 'Position', [.8314 .11 .03 .8])
     
    
end

%% Plotting eigen basis 
% 
% for i=1:size(eigVec,2)
%     contourDomainPlot(fem,idPart,eigVec(:,i),1);
%     print(sprintf('Basis%i',i),'-dpng') 
% end
