%%Partial measurement data case implementaion for hinge--before the
%%function implementation
close all;
clear ;
dbstop if error;

run 'C:\Users\babu_m\Documents\GitHub\source\initFiles.m';
eps=10E-16;
fem=femInit();
fileName{1}=strcat(inpFiles,'\Remesh_hinge.inp');%\top_hat.inp');%\halo.inp');%
% fem=importMesh(fem,fileName{1});%'top_hat.inp'); %use only for top hat
fem=importMultiMesh(fem, fileName);%has provisions for domain separation
% within code
fem=femPreProcessing(fem);
domainID=1;
idPart=domainID;
coordi=fem.xMesh.Node.Coordinate;
nodeIdDomain=fem.Domain(idPart).Node;
nodeCoord=fem.xMesh.Node.Coordinate(nodeIdDomain,:); % coordinates for the domain
%load pre-defined part regions
load('hingeRegionsCoarse.mat');% contains variable nodeIDoutRect
%load optimised error process covariance parameters
load('optimisedHypParm.mat'); %contains struct hyp in gp function format(i.e. hyp.cov=log(parameters)) for matern covariance

%% for coarse mesh
devPatterns= load('hingeDevArSimLocCoarse.mat');%  %hingeDevArSimLoc
devPatterns=devPatterns.femDevDomain;
% % % testing the effect of not centering
devCenterd1=devPatterns;
% devCenterd1=bsxfun(@minus,devPatterns,mean(devPatterns)); % to calculate the rmse and for modelling deviation from mean
% devSelected= devPatterns(:,iMnp);
devCenterd=devCenterd1;%devCenterd1+randn(size(devPatterns))*10e-3;%devPatterns;x=z*sigma+mu;; % set equal to input in a coarse mesh

% to obtain major deviation patterns
[U ,singVal,eigVec]=svds(devCenterd,15);
eigVal=diag(singVal).^(2)./size(devCenterd,2);

pc=devCenterd*eigVec;%plot(pc);

x=pc*eigVec';
diff=devCenterd-x;


%% the kalman filter notations are from Cressie and wikle (biometrica)
%Z_t=eigVec*a_t+\nu+\eps
%a_t=(eigVec'*B)*a_(t-1)+eigVec'*\eta
%a_t=H*a_(t-1)+eigVec'*\eta

% estimation of covariances c_z
[n, t]=size(devCenterd);
% zero lag covariance
gMean=sum(sum(devCenterd))/(n*t);
C0z=((devCenterd-gMean)'*(devCenterd-gMean))./t;
% lag one covariance
devCenterdT=devCenterd(2:end,:);
devCenterdTminus1=devCenterd(1:end-1,:);
C1z=((devCenterdT-gMean)'*(devCenterdTminus1-gMean))./(t-1);
C1z = nearestSPD(C1z);
% measurement error covariance
sigmaMes=1E-4;                               % measurement error variance
R=sigmaMes.*eye(size(C1z));
% coloured error process covariance matrix
covFunc = {@covMaternard,3};   % Matern class d=3
V=eye(size(R))*1e-4;%feval(covFunc{:},hyp.cov,nodeCoord); % optimised values obtained from average of maximum likelihood values of all samples
% V=feval(@covMaternard,3,hyp.cov,nodeCoord); % this method to evaulate cov func also works
% zero lag y process
C0y=C0z-R;
C0y = nearestSPD(C0y);
% B matrix ie. the matrix that specifies autocorrelation b/w a_t and a_(t-1)
J=eigVec'; % I take eigVec'eigVec to be identity matrix as eigVec is orthogonal (eigVec corresponds to phi the eigen vector matrix)
B=(C1z*J')/(J*(C0z-R-V)*J');
% estimation of JQJ' matrix ie. the state error covariance matrix
H=eigVec'*B;
covErr=J*(C0z-R-V)*J'-J*C1z*J'*H'-H*J*C1z'*J'+H*J*(C0z-R-V)*J'*H';
covErr = nearestSPD(covErr);


%% kalman recursion

zt=devCenterd'; % measurements
tol=2;    % is a trickey values affects the adaptivity a lot..currently works well when tol is set to 2
countPart=1;
maxSnap=4;
numberCompleteMeasure= 15; % to stabalize the prediction
subPlotRos=2;
subPlotCol=2;
cmin=-3;
cmax=3;

% at time zero
a00=pc(1,:)';
p00=eye(size(eigVec,2));
ytp1=zeros(size(nodeCoord,1),size(devCenterd,1));
ytp1=ytp1';
ytt=zeros(size(ytp1));
VarYtt=zeros(size(ytp1));
VarYtp1=zeros(size(ytp1));


% to stabalise the filtering equations
for i=1:numberCompleteMeasure
    
    if i==1
        attm1=a00;
        pttm1=p00;
    else
        %prediction
        attm1=H*att;
        pttm1=H*ptt*H'+covErr;
    end
    
    %observation
    kt=(pttm1*eigVec')/(R+V+eigVec*pttm1*eigVec');
    att=attm1+kt*(zt(:,i)-eigVec*attm1);
    ptt=pttm1-kt*eigVec*pttm1;
    
    ytt(i,:)=att'*eigVec';   % taking transpose to maintain convention of deviation matrix(i.e. rep x nodes)
    VarYtt(i,:)=diag(eigVec*ptt*eigVec'+(V+R));
    ytp1(i,:)=att'*H'*eigVec';  % taking transpose to maintain convention of deviation matrix
    VarYtp1(i,:)=diag(eigVec*(H*ptt*H'+covErr)*eigVec'+(V+R));
    
end
%%

for i=numberCompleteMeasure+1:size(zt,2) % i for each part
    
    nSnaps=1;
    allMesReg=[];
    nRegions=1:size(nodeIDoutRect,2);
    fig1=figure;
    fig2=figure;
    fig3=figure;
    rmseT=zeros(maxSnap,1);
    
    while ~isempty(nRegions)&&(nSnaps<=maxSnap)
        selIndex=zeros(length(nRegions),1);
        
        % prediction based on t-1 to predict first region to measure
        for j=1:length(nRegions)
            if nSnaps==1  % when no region is measured use previous prediction
                regionID=nRegions(j);
                regionIndex=find(nodeIDoutRect(:,regionID));%% finds selected nodes
                devRegion=ytp1(i-1,regionIndex);
                varRegion=VarYtp1(i-1,regionIndex);
                % calculating the probability of being out of tolerance
                probOut=1-cdf('Normal',tol,devRegion,sqrt(varRegion))+cdf('Normal',-tol,devRegion,sqrt(varRegion));
                probOut(probOut<=eps)=1; % to avoid the infinity case and take the limit
                selIndex(j)=sum(-probOut.*log(probOut))/length(regionIndex);  % remember to normalise to counter the effect of number of points in the region
            else % when some regions are measured use updated predictionregionID=nRegions(j);
                regionID=nRegions(j);
                regionIndex=find(nodeIDoutRect(:,regionID));%% finds selected nodes
                devRegion=ytt(i,regionIndex);
                varRegion=VarYtt(i,regionIndex);
                % calculating the probability of being out of tolerance
                probOut=1-cdf('Normal',tol,devRegion,sqrt(varRegion))+cdf('Normal',-tol,devRegion,sqrt(varRegion));
                probOut(probOut<=eps)=1;
                selIndex(j)=sum(-probOut.*log(probOut))/length(regionIndex);  % remember to normalise to counter the effect of number of points in the region
                
            end
        end
        
        
        %mesurement vector concatenation
        if nSnaps==1
            [~,indexMax]=max(selIndex);
            selRegion=nRegions(indexMax);
            allMesReg=[allMesReg;selRegion]; % all the measured regions
            tempV=rearrangeCovarianceMat(V,nodeIDoutRect(:,selRegion));
            tempR=rearrangeCovarianceMat(R,nodeIDoutRect(:,selRegion));
            regionIndex=find(nodeIDoutRect(:,selRegion));%% finds selected nodes
            tempEV=eigVec(regionIndex,:);
            tempZ= zt(regionIndex,i);
            
            tempKt=(pttm1*tempEV')/(tempV+tempR+tempEV*pttm1*tempEV');
            att=attm1+tempKt*(tempZ-tempEV*attm1);
            ptt=pttm1-tempKt*tempEV*pttm1;
            
            ytt(i,:)=att'*eigVec';   % taking transpose to maintain convention of deviation matrix(i.e. rep x nodes)
            VarYtt(i,:)=diag(eigVec*ptt*eigVec'+(V+R));
            ytp1(i,:)=att'*H'*eigVec';  % taking transpose to maintain convention of deviation matrix
            VarYtp1(i,:)=diag(eigVec*(H*ptt*H'+covErr)*eigVec'+(V+R));
            
            nRegions(indexMax)=[];
            nSnaps=nSnaps+1;
            
        else
            [~,indexMax]=max(selIndex);
            selRegion=nRegions(indexMax);
            allMesReg=[allMesReg;selRegion]; % all the measured regions
            allPoints=zeros(length(nodeIdDomain),1);
            
            for j=1:length(allMesReg)
                allPoints=allPoints|nodeIDoutRect(:,allMesReg(j));
            end
            
            tempV=rearrangeCovarianceMat(V,allPoints);
            tempR=rearrangeCovarianceMat(R,allPoints);
            regionIndex=find(allPoints);%% finds selected nodes
            tempEV=eigVec(regionIndex,:);
            tempZ= zt(regionIndex,i);
            
            tempKt=(pttm1*tempEV')/(tempV+tempR+tempEV*pttm1*tempEV');
            att=attm1+tempKt*(tempZ-tempEV*attm1);
            ptt=pttm1-tempKt*tempEV*pttm1;
            
            ytt(i,:)=att'*eigVec';   % taking transpose to maintain convention of deviation matrix(i.e. rep x nodes)
            VarYtt(i,:)=diag(eigVec*ptt*eigVec'+(V+R));
            ytp1(i,:)=att'*H'*eigVec';  % taking transpose to maintain convention of deviation matrix
            VarYtp1(i,:)=diag(eigVec*(H*ptt*H'+covErr)*eigVec'+(V+R));
            
            nRegions(indexMax)=[];
            nSnaps=nSnaps+1;
            
        end
        
        % plot for prediction after each measurement
        set(0, 'currentfigure', fig1)
        plotAxis1=subplot(subPlotRos,subPlotCol,nSnaps-1);
        contourDomainPlot(fem,idPart,ytt(i,:),0,plotAxis1)
        colormap('jet')
        hold on
        points2plot=nodeCoord(regionIndex,:);
        scatter3(points2plot(:,1),points2plot(:,2),points2plot(:,3),'+k');
        hold off
        caxis([cmin,cmax]);
        set(0, 'currentfigure', fig2)
        plotAxis2=subplot(subPlotRos,subPlotCol,nSnaps-1);
        contourDomainPlot(fem,idPart,VarYtt(i,:),0,plotAxis2)
        colormap('jet')
        
        %plot without measured region
        set(0, 'currentfigure', fig3)
        plotAxis1=subplot(subPlotRos,subPlotCol,nSnaps-1);
        contourDomainPlot(fem,idPart,ytt(i,:),0,plotAxis1)
        caxis([cmin,cmax]);
        colormap('jet')
        % RMSE calculation for the step
        rmseTnoisy(nSnaps)=sqrt(sum((bsxfun(@minus,devCenterd(i,:),ytt(i,:))).^2,2)./size(devCenterd,2)); % rmse with noisy input
%         rmseT(nSnaps)=sqrt(sum((bsxfun(@minus,devCenterd1(i,:),ytt(i,:))).^2,2)./size(devCenterd1,2));  % rmse without noisy input
        
        
    end
    
    % Updating the filtering equations for next time step
    attm1=H*att;
    pttm1=H*ptt*H'+covErr;
    
    %% plotting 
    
%     %plotting RMSE
%     h=figure;
%     plot(rmseT(2:end));
%     stringTitle=sprintf('Noise free RMSE for prediction in time t after each measurements');
%     ax=h.Children;
%     ax.XTick=1:maxSnap;
%     title(stringTitle);
%     xlabel('Number of measurements');
%     ylabel('RMSE in mm');
    
    %plotting RMSE
    h=figure;
    plot(rmseTnoisy(2:end));
    ax=h.Children;
    ax.XTick=1:maxSnap;
    stringTitle=sprintf('RMSE for prediction in time t after each measurements');
    title(stringTitle);
    xlabel('Number of measurements');
    ylabel('RMSE in mm');
    
    %plotting original Deviation
    contourDomainPlot(fem,idPart,zt(:,i),1)
    stringTitle=sprintf('Original variation pattern');
    title(stringTitle);
    caxis([cmin,cmax]); 
    colormap('jet')
    
    %plotting prediction for next time step
    contourDomainPlot(fem,idPart,ytp1(i,:),1)
    stringTitle=sprintf('Prediction for next time step');
    title(stringTitle);
    caxis([cmin,cmax]); 
    colormap('jet')
    
    % setting title to figures
    set(0, 'currentfigure', fig1)
    stringTitle=sprintf('Sequence of adaptive measurements based on information gain');
    figtitle(stringTitle);
    
    set(0, 'currentfigure', fig2)
    % to set single colourbar
    colorbar
    h=colorbar;
    set(get(h,'Ylabel'),'string','mm^{2}','fontweight','bold')
    set(h, 'Position', [.8314 .11 .03 .8150])
    for kk=1:maxSnap
        ax=subplot(subPlotRos,subPlotCol,kk);
        pos=get(ax, 'Position');
        set(ax, 'Position', [pos(1) pos(2) 0.85*pos(3) pos(4)]);
    end
    % setting title before colourbar messes the title setting function
    stringTitle=sprintf('Variance of predicted measurements after each measurement');
    figtitle(stringTitle);
    
    set(0, 'currentfigure', fig3)  
    % to set single colour bar
    colorbar
    h=colorbar;
    set(get(h,'Ylabel'),'string','mm','fontweight','bold')
    set(h, 'Position', [.8314 .11 .03 .8150])
    for kk=1:maxSnap
        ax=subplot(subPlotRos,subPlotCol,kk);
        pos=get(ax, 'Position');
        set(ax, 'Position', [pos(1) pos(2) 0.85*pos(3) pos(4)]);
    end
    stringTitle=sprintf('Succesive predictions for entire part after each measurement');
    figtitle(stringTitle);
end
