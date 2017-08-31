%% basis function generation (EOF or principal components)and kalman filter implementaiton
% has maxlikelihood estimation for the covariance of the residuals of
% the kalman filter equations, does not take \nu into account 
close all; clear ;

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


%% eigen function calculation

% % for dense mesh
% [ Pp, Mnp,iMnp  ] = uniformKeyPointSelection( fem,idPart,20 );
% % [C,iMnp,ic] = unique(Mnp,'rows'); % has not been implemented yet
% devPatterns= load('hingeDevArSimLoc.mat');%
% devPatterns=devPatterns.femDevDomain;
% devCenterd1=bsxfun(@minus,devPatterns,mean(devPatterns)); % to calculate the rmse, used nowhere else
% devSelected= devPatterns(:,iMnp);
% devCenterd=bsxfun(@minus,devSelected,mean(devSelected));


%% for coarse mesh
% [ Pp, Mnp,iMnp  ] = uniformKeyPointSelection( fem,idPart,20 );
devPatterns= load('hingeDevArSimLocCoarse.mat');%  %hingeDevArSimLoc
devPatterns=devPatterns.femDevDomain;
devCenterd1=bsxfun(@minus,devPatterns,mean(devPatterns)); % to calculate the rmse, used nowhere else
% devSelected= devPatterns(:,iMnp);
devCenterd=devCenterd1; % set equal to input in a coarse mesh


% percentVar=98/100;
% [U,eigVal,V]=getEigenTurnc(devCenterd,percentVar);

[U ,singVal,V]=svds(devCenterd,15);
eigVal=diag(singVal).^(2)./size(devCenterd,2);

pc=devCenterd*V;%plot(pc);

x=pc*V';
diff=devCenterd-x;

%% fitting the covariance parameters for the error process
trainX=nodeCoord;
covParam=zeros(size(diff,1),size(trainX,2)+1);% for the matern and most other covariance function
likParam=zeros(size(diff,1),1);

%% max likelihood calculation
% for i=1:size(diff,1)
%     trainY=diff(i,:)';
%     L = rand(size(trainX,2),1)*400;%[10;10;30];%
%     Liso=100;
%     m0 = {'meanZero'};  hyp.mean = [];      % no hyperparameters are needed
%     % m0 = {@meanConst};  hyp.mean = 1;  % also function handles are possible
%     sn=0.000001; % noise standard dev
%     likfunc = {'likGauss'};   hyp.lik = log(sn);  % gaussian likelihood
%     % mu = 1.0; s2 = 0.01^2;
%     % pg = {@priorGauss,mu,s2};                          % Gaussian prior
%     sf=0.01;           % latent function standard dev
%     % covFunc = {@covSEard};   hyp.cov = log([L;sf]);
%     % covFunc={@covSEiso};  hyp.cov = log([Liso;sf]);    % isotropic Gaussian
%     covFunc = {@covMaternard,3};  hyp.cov = log([L;sf]); % Matern class d=3
%     % covFunc = {'covMaterniso',3}; hyp.cov = log([Liso;sf]); % Matern class d=3
%     hypOpt= minimize(hyp, @gp, -50,  @infExact, m0, covFunc, likfunc, trainX,trainY );
%     covParam(i,:)=hypOpt.cov;
%     likParam(i,:)=hypOpt.lik;
%     
% end
% covParam=exp(covParam);
% likParam=exp(likParam);
% 
% % using the mean values as the optimised value
% covMean=mean(covParam);
% likMean=mean(likParam);
% % simple krigging the residuals
% hyp.cov=log(covMean);
% hyp.lik=log(likMean);
% % [mpredicted, s2] = gp(hyp, @infExact, m0, covFunc, likfunc, xt, yt, z);


figure
imagesc(diff);
figure;
rmsePC=sqrt(sum((diff).^2,2)./size(devCenterd1,2));
% the error before and after kalman filter remains same because the
% correlated error process (\nu   ) is not taken care while deriving the kalman
% filter equations
plot(rmsePC);

%% the kalman filter notations are from Cressie and wikle (biometrica)

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
% zero lag y process
C0y=C0z-R;
C0y = nearestSPD(C0y);
% B matrix ie. the matrix that specifies autocorrelation b/w a_t and a_(t-1)
J=V'; % I take V'V to be identity matrix as v is orthogonal (V corresponds to phi the eigen vector matrix)
B=(C1z*J')/(J*(C0z-R)*J');
% estimation of JQJ' matrix ie. the error covariance matrix
H=V'*B;
covErr=J*(C0z-R)*J'-J*C1z*J'*H'-H*J*C1z'*J'+H*J*(C0z-R)*J'*H';
covErr = nearestSPD(covErr);

eigVec=V;
%% kalman recursion

% at time zero
a00=pc(1,:)';
p00=eye(size(V,2));
ytp1=zeros(size(nodeCoord,1),size(devCenterd,1));
ytp1=ytp1';
ytp=zeros(size(ytp1));
zt=devCenterd';

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
    
    
    %observation
    kt=(pttm1*V')/(R+V*pttm1*V');
    att=attm1+kt*(zt(:,i)-V*attm1);
    ptt=pttm1-kt*V*pttm1;
    
    ytp(i,:)=att'*eigVec';   % taking transpose to maintain convention of deviation matrix(i.e. rep x nodes)
    ytp1(i,:)=att'*H'*eigVec';  % taking transpose to maintain convention of deviation matrix
    
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
    plotAxis1=subplot(2,2,i);
    contourDomainPlot(fem,idPart,devCenterd1(i,:),plotAxis1)
    %     caxis([-2 2]);
    
    set(0, 'currentfigure', fig2)
    plotAxis2=subplot(2,2,i);
    contourDomainPlot(fem,idPart,ytp(i,:),plotAxis2)
    %     caxis([-2 2]);
    
    set(0, 'currentfigure', fig3)
    plotAxis3=subplot(2,2,i);
    contourDomainPlot(fem,idPart,ytp1(i,:),plotAxis3)
    %     caxis([-1.5 1.5]);
    
    
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
rmseT=sqrt(sum((bsxfun(@minus,devCenterd1(2:end,:),ytp1(1:end-1,:))).^2,2)./size(devCenterd1,2));
plot(rmseT);
stringTitle=sprintf('RMSE for prediction for time (t+1) after measuring %d points in time t',length (nodeIdDomain));
title(stringTitle);
xlabel('instance');
ylabel('RMSE in mmm');