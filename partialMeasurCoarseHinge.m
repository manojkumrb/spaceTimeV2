%%Partial measurement data case implementaion for hinge
close all;
clear ;

run 'C:\Users\babu_m\code\source\initFiles.m';

fem=femInit(source);
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
% devCenterd1=devPatterns;
devCenterd1=bsxfun(@minus,devPatterns,mean(devPatterns)); % to calculate the rmse and for modelling deviation from mean
% devSelected= devPatterns(:,iMnp);
devCenterd=devCenterd1; % set equal to input in a coarse mesh

% to obtain major deviation patterns
[U ,singVal,EV]=svds(devCenterd,15);
eigVal=diag(singVal).^(2)./size(devCenterd,2);

pc=devCenterd*EV;%plot(pc);

x=pc*EV';
diff=devCenterd-x;


%% the kalman filter notations are from Cressie and wikle (biometrica)
%Z_t=EV*a_t+\nu+\eps
%a_t=(EV'*B)*a_(t-1)+EV'*\eta
%a_t=H*a_(t-1)+EV'*\eta

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
sigmaMes=1E-12;                               % measurement error variance
R=sigmaMes.*eye(size(C1z));
% coloured error process covariance matrix
covFunc = {@covMaternard,3};   % Matern class d=3
V=eye(size(R))*1e-10;%feval(covFunc{:},hyp.cov,nodeCoord); % optimised values obtained from average of maximum likelihood values of all samples
% V=feval(@covMaternard,3,hyp.cov,nodeCoord); % this also works
% zero lag y process
C0y=C0z-R;
C0y = nearestSPD(C0y);
% B matrix ie. the matrix that specifies autocorrelation b/w a_t and a_(t-1)
J=EV'; % I take EV'EV to be identity matrix as EV is orthogonal (EV corresponds to phi the eigen vector matrix)
B=(C1z*J')/(J*(C0z-R-V)*J');
% estimation of JQJ' matrix ie. the state error covariance matrix
H=EV'*B;
covErr=J*(C0z-R-V)*J'-J*C1z*J'*H'-H*J*C1z'*J'+H*J*(C0z-R-V)*J'*H';
covErr = nearestSPD(covErr);
eigVec=EV;

%% kalman recursion

% at time zero
a00=pc(1,:)';
p00=eye(size(EV,2));
ytp1=zeros(size(nodeCoord,1),size(devCenterd,1));
ytp1=ytp1';
ytp=zeros(size(ytp1));
zt=devCenterd';



%%
tol=1;
nRegions=1:size(nodeIDoutRect,2);
count=1;
limit=4;
%%
for i=1:size(zt,2)
    
    %% kalman filter
    if i==1
        attm1=a00;
        pttm1=p00;
    else
        %prediction
        attm1=H*att;
        pttm1=H*ptt*H'+covErr;
    end
    
    if i==9
        
        while ~isempty(nRegions)&&(count<=limit)
            selIndex=zeros(length(nRegions),1);
            
            % prediction based on t-1 to predict first region to measure
            for j=1:length(nRegions)
                
                tempV=rearrangeCovarianceMat(V,nodeIDoutRect(:,j));
                tempR=rearrangeCovarianceMat(R,nodeIDoutRect(:,j));
                regionIndex=find(nodeIDoutRect(:,j));%% finds selected nodes
                tempEV=EV(regionIndex,:);
                tempZ= zt(regionIndex,i);               
                
                tempKt=(pttm1*tempEV')/(tempR+tempV+tempEV*pttm1*tempEV');
                att=attm1+tempKt*(tempZ-tempEV*attm1);
                ptt=pttm1-kt*EV*pttm1;
                
                ytp(i,:)=att'*eigVec';   % taking transpose to maintain convention of deviation matrix(i.e. rep x nodes)
                VarYtt=diag(eigVec*ptt*eigVec'+(V+R));
                ytp1(i,:)=att'*H'*eigVec';  % taking transpose to maintain convention of deviation matrix
                VarYtp1=diag(eigVec*(H*ptt*H'+covErr)*eigVec'+(V+R));
                % calculating the probability of being out of tolerance
                probOut=1-cdf('Normal',tol,ytp(i,:)',sqrt(VarYtt))+cdf('Normal',-tol,ytp(i,:)',sqrt(VarYtt));
                selIndex(j)=sum(-probOut.*log(probOut));
                
            end
            [~,indexMax]=max(selIndex);
            tempV=rearrangeCovarianceMat(V,nodeIDoutRect(:,indexMax));
            tempR=rearrangeCovarianceMat(r,nodeIDoutRect(:,indexMax));                

            tempKt=(pttm1*tempEV')/(tempR+tempR+tempEV*pttm1*tempEV');
            att=attm1+tempKt*(tempZ-tempEV*attm1);
            ptt=pttm1-kt*EV*pttm1;

            ytp(i,:)=att'*eigVec';   % taking transpose to maintain convention of deviation matrix(i.e. rep x nodes)
            VarYtt=diag(eigVec*ptt*eigVec'+(V+R));
            ytp1(i,:)=att'*H'*eigVec';  % taking transpose to maintain convention of deviation matrix
            VarYtp1=diag(eigVec*(H*ptt*H'+covErr)*eigVec'+(V+R));
            % calculating the probability of being out of tolerance
            probOut=1-cdf('Normal',tol,ytp(i,:)',sqrt(VarYtp1))+cdf('Normal',-tol,ytp(i,:)',sqrt(VarYtp1));
            selIndex(j)=sum(-probOut.*log(probOut));
                      
            
            
            
            nRegions(indexMax)=[];
            count=count+1;
            clear score
            
        end
        
    end
    %observation
    kt=(pttm1*EV')/(V+R+EV*pttm1*EV');
    att=attm1+kt*(zt(:,i)-EV*attm1);
    ptt=pttm1-kt*EV*pttm1;
    
    ytp(i,:)=att'*eigVec';   % taking transpose to maintain convention of deviation matrix(i.e. rep x nodes)
    VarYtt=diag(eigVec*ptt*eigVec'+(V+R));
    ytp1(i,:)=att'*H'*eigVec';  % taking transpose to maintain convention of deviation matrix
    VarYtp1=diag(eigVec*(H*ptt*H'+covErr)*eigVec'+(V+R));
end
%
% figure('comparision of predicted and actual values')
% for i=1:size(devCenterd,1)-1
%     plot(ytp1(:,i)'-zt(i+1,:));
%     hold on
% end

fig1=figure;

fig2=figure;

fig3=figure;


for i=1:4
    set(0, 'currentfigure', fig1)
    %         p=i+10; % to plot values other than the initila deviation values
    
    plotAxis1=subplot(2,2,i);
    contourDomainPlot(fem,idPart,devCenterd1(p,:),plotAxis1)
    %     caxis([-2 2]);
    
    set(0, 'currentfigure', fig2)
    plotAxis2=subplot(2,2,i);
    contourDomainPlot(fem,idPart,ytp(p,:),plotAxis2)
    %     caxis([-2 2]);
    
    set(0, 'currentfigure', fig3)
    plotAxis3=subplot(2,2,i);
    contourDomainPlot(fem,idPart,ytp1(p,:),plotAxis3)
    %     caxis([-1.5 1.5]);
    %%% simple krigging the residuals
    %     hyp.cov=log(covMean);
    %     hyp.lik=log(likMean);
    %     [mpredicted, s2] = gp(hyp, @infExact, m0, covFunc, likfunc, xt, yt, z);
end

set(0, 'currentfigure', fig1);
figtitle('Auto-correlated deviations generated through simulation','fontweight','bold');
set(0, 'currentfigure', fig2);
figtitle('Predicted deviation for the entire part at time t','fontweight','bold');
set(0, 'currentfigure', fig3)
figtitle('Predicted deviation for the entire part at time (t+1)','fontweight','bold');

figure;
rmseT=sqrt(sum((bsxfun(@minus,devCenterd1,ytp)).^2,2)./size(devCenterd1,2));
plot(rmseT);
stringTitle=sprintf('RMSE for prediction in time t after measuring %d points',length (nodeIdDomain));
title(stringTitle);
xlabel('instance');
ylabel('RMSE in mmm');

figure;
rmseTp1=sqrt(sum((bsxfun(@minus,devCenterd1(2:end,:),ytp1(1:end-1,:))).^2,2)./size(devCenterd1,2));
plot(rmseTp1);
stringTitle=sprintf('RMSE for prediction for time (t+1) after measuring %d points in time t',length (nodeIdDomain));
title(stringTitle);
xlabel('instance');
ylabel('RMSE in mmm');

