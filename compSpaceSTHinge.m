% code for comparision of space time KF and krigging


close all;clear all;
run 'C:\Users\babu_m\Documents\GitHub\source\initFiles.m';
path(path,'\\megara\dlmr\Manoj\all data sensitivity inner'); % all data files in megara
load('hingeDevations.mat'); % contains variable dev
load('hingeDevArSimLocCoarse.mat'); % contains variable femDevDomain
load('seqPedStructHingeCoarseWithRmseTp1_25EigWithKrig.mat'); % contains structure seqPred %

devSmooth       =femDevDomain;
% load('hatLowerDev.mat'); % contains variable devHatLower

% load pre-defined part regions
partRegions     =load('hingeRegionsCoarse.mat');% contains variable nodeIDoutRect
nodeIDoutRect   =partRegions.nodeIDoutRect;


fem             =femInit();
fileName{1}     =strcat(inpFiles,'\Remesh_hinge.inp');%\hinge.inp');%\top_hat.inp');%\halo.inp');%
fem             =importMultiMesh(fem,fileName);%'top_hat.inp'); % has provisions for domain separation within code
totalNodes      =size(fem.xMesh.Node.Coordinate,1);
fem             =femPreProcessing(fem);
domainID        =1;
idPart          =domainID;

nodeIdDomain    =fem.Domain(idPart).Node;
nodeCoord       =fem.xMesh.Node.Coordinate(nodeIdDomain,:); % coordinates for the domain
nodeCoord2      =nodeCoord(:,[1,3]); % make node 2d

%% learning parameters
%%%%%%%%%%%%%%%%%%%%%%%%%
sourcePattern   =18;%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
x               = nodeCoord; % create a data set
L               = rand(size(x,2),1)*10;% [50;50;250];%

%% mean and covariance setting
m0              ={'meanZero'};  hyp.mean = [];                      % no hyperparameters are needed
sn              =0.001;                                             % noise standard dev
likfunc         ={'likGauss'};   hyp.lik = log(sn);                 % gaussian likelihood
sf              =1;                                                 % latent function standard dev
lScale          =[60;50;100]*rand(1);                               %rand(D,1);% characteristic length scale
covMatern       ={@covMaternard,5};  hypMat=log([lScale;sf]);      % Matern class d=5
% covMatern       ={@covSEard};  hypMat=log([lScale;sf]);      % modified matern to sq exp
csu             ={'covSum',{covMatern}}; hypsu = hypMat;         % noise linear, quadratic, periodic
% cra = {'covRQard'}; al = 2; hypra = log([L;sf;al]); % ration. quad.
% csu = {'covSum',{cra}}; hypsu = hypra; % noise linear, quadratic, periodic
hyp.cov         =hypsu;
nodeLimit       =200;
xLimited        =sampleNodesDomain(fem,domainID,nodeLimit);
covFuncI        ={@covFITC, csu,xLimited};

%% log marginal likelihood optimization
% hyp = minimize(hyp, @gp, -50,  @infFITC, m0, covFuncI, likfunc,  x, devSmooth(sourcePattern,:)); % for covFITC
hyp             = minimize(hyp, @gp, -30, @infExact, m0, csu, likfunc, x, devSmooth(sourcePattern,:)'); % for regular covariance calculation

%%
numberCompleteMeasure   =10;
maxSnap                 =length(seqPred(end).Snap);
nPattern                =size(devSmooth,1);
gprPred(nPattern).Snap(1).Dev=[];


for i=numberCompleteMeasure+1:nPattern % i for each part
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    inputPattern=i;%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    rmse=zeros(maxSnap,1);
    
    for j=1:maxSnap
        
        %% Selecting region with measurement
        allMesReg       =seqPred(inputPattern).Snap(j).Region;
        allPoints       =zeros(size(nodeIDoutRect,1),1);
        for k =1:length(allMesReg)
            allPoints   =allPoints|nodeIDoutRect(:,allMesReg(k));
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
        [mpredicted, s2] = gp(hypPred, @infExact, m0, csu, likfunc, xt, yt', z);
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
cmin        =-3;
cmax        =3;
nInstances  =5; % number of parts to plot
cDelta=0.25;     % step size for contourmap


for i=numberCompleteMeasure+2:numberCompleteMeasure+2+nInstances
    
    
    rmseST=seqPred(i).Snap(maxSnap).RmseT;
    rmseSTp1=seqPred(i-1).Snap(maxSnap).rmseTp1;
    rmseST=[rmseSTp1;rmseST];
    rmseS=gprPred(i).Snap(maxSnap).Rmse;
    %     rmseS=[0;rmseS];
    
%     
%     for kk=1:maxSnap-3
%         
%         devST=seqPred(i).Snap(kk).Dev;
%         devS=gprPred(i).Snap(kk).Dev;
%         
%         var=seqPred(i).Snap(kk).Var;
%         allMesReg=seqPred(i).Snap(kk).Region;
%         
%         fig1        =figure;
%         % plot for prediction after each measurement
%         set(0, 'currentfigure', fig1)
%         plotAxis1   =subplot(subPlotRos,subPlotCol,1);
%         contourDomainPlot(fem,idPart,devSmooth(i,:),0,plotAxis1)
%         
%         plotAxis1   =subplot(subPlotRos,subPlotCol,2);
%         contourDomainPlot(fem,idPart,devS,0,plotAxis1)
%         
%         plotAxis1   =subplot(subPlotRos,subPlotCol,3);
%         contourDomainPlot(fem,idPart,devST,0,plotAxis1)
%         
%         contourcmap('parula',cDelta);
%         %         fig1=makeSubplotCompact(fig1,subPlotRos, subPlotCol);
%         
%     end
    
    %%% Plot RMSE
    h=figure;
    myLine2=plot(rmseST);
    hold all;
    myLine=plot(2:length(rmseST),rmseS);
    
    
    
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
    
    
end

% % average Rmse diff
% 
% 
% for kk=1:maxSnap
%     diff=0;
%     rmseST=0; 
%     rmseS=0;
%     
%     for i=numberCompleteMeasure+1:nPattern
%         
%         rmseST=rmseST+seqPred(i).Snap(kk).RmseT(kk);
%         rmseS=rmseS+gprPred(i).Snap(kk).Rmse(kk);
%         
%     end
%     
%     delta(kk)=((rmseS-rmseST)/rmseS)*100;
%     
% end
% fig=figure;
% hh=plot(delta(1:6));
% fig=changeAxesLooks(fig,{'Average improvement in RMSE with Spatio-temporal prediction','compared to spatial prediction '},...
%     'Number of measurements','Percentage improvement');
% hh.LineWidth=1.5;

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

totalInstances= nPattern-numberCompleteMeasure+1;

delta=((rmseS-rmseST)./rmseS)*100;
fig=figure;
hh=plot(delta);
fig=changeAxesLooks(fig,{'Average improvement in RMSE with Spatio-temporal prediction','compared to spatial prediction '},...
    'Number of measurements','Percentage improvement');

fig=figure;
hh=plot(1:length(rmseS),rmseS./totalInstances,1:length(rmseST),rmseST./totalInstances);
fig=changeAxesLooks(fig,{'Average RMSE of Spatio-temporal prediction','compared to spatial prediction '},...
    'Number of measurements','RMSE in mm');

set(hh,'linewidth',1.5);

legendflex(hh, {'Spatial prediction','Spatio-temporal prediction'}, ...
        'anchor', {'ne','ne'}, ...
        'buffer', [0 0], ...
        'nrow',   2,...
        'fontsize',10, ...
        'xscale', 0.66, ...
        'box','on');