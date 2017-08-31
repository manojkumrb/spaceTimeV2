%% basis function generation (EOF or principal components)
% made changes to include edge nodes
close all; clear ;

run 'C:\Users\babu_m\Documents\GitHub\source\initFiles.m';

fem=femInit();
fileName{1}=strcat(inpFiles,'\door_inner.inp');%\hinge.inp');%\3dTriaHinge.inp');%\halo.inp');%
% fem=importMesh(fem,fileName{1});%'top_hat.inp'); % has provisions for domain separation within code
fem=importMultiMesh(fem, fileName);
fem=femPreProcessing(fem);
totalNodes=size(fem.xMesh.Node.Coordinate,1);
domainID=1;
idPart=domainID;
coordi=fem.xMesh.Node.Coordinate;
nodeIdDomain=fem.Domain(idPart).Node;
nodeCoord=fem.xMesh.Node.Coordinate(nodeIdDomain,:); % coordinates for the domain

%% %% graphical selction
[ ~, ~,iMnp  ] = uniformKeyPointSelection( fem,idPart,45);
[boundNodeID]=getBoundaryNodes(fem,domainID, 150); % outer
[boundNodeID2]=getBoundaryNodes(fem,domainID, 40);  % halo
% [boundNodeID3]=getBoundaryNodes(fem,domainID, 40);  % below
manualSelect=[5230;24031;5030;4698;4862;4374;4531;4729;4201;4340;...
	4511;4039;4171;4337;3838;3982;4183;3677;3756;4008;3646;3583;904;837];
% manualSelect2=[52739;52831;52875;53030;53123;69790;75566;75875;54576;63270;...
% 	63317;63355;54377;60645;49175;64910;64766;65396;65584;65706;61196;61340];
iMnp=[iMnp; boundNodeID;boundNodeID2;manualSelect];
uniqeMnp= unique(iMnp);
Mnp=nodeCoord(uniqeMnp,:);
selNodes.ID=uniqeMnp;
selNodes.Coord=Mnp;

%% checking selected nodes
meshplotAxisDefined(fem,1);
hold all
scatter3(Mnp(:,1),Mnp(:,2),Mnp(:,3),'*');

%% loading from predefined selection
selNodes=load('doorInnerSelNodes3.mat');
selNodes=selNodes.selNodes;
Mnp=selNodes.Coord;
iMnp=selNodes.ID;
% iMnp=dsearchn(nodeCoord,Mnp);


%%
devPatterns= load('devAutoCorrInner_dense.mat');%flatPlateDevArLocators.mat');%
devPatterns=devPatterns.femDevDomain;%femDevArFlatPlate;
devCenterd1=bsxfun(@minus,devPatterns,mean(devPatterns)); % to calculate the rmse, used nowhere else
devSelected= devPatterns(:,iMnp);
devCenterd=bsxfun(@minus,devSelected,mean(devSelected));

% percentVar=99/100;
% [U,eigVal,V]=getEigenTurnc(dev_centerd,percentVar);

% [U ,singVal,V]=svds(devCenterd,20);
% eigVal=diag(singVal).^(2)./size(devCenterd,2);

% 
V=load('eigVecInner3_400.mat');
V=V.V;
pc=devCenterd*V;%plot(pc);
x=pc*V';
diff=devCenterd-x;
figure
imagesc(diff);


%% the kalman filter notations are from Cressie and wikle (biomtrica)

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
sigmaMes=1E-12;                               % measurement variance
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


%% interpolation of eigen vectors

% % removing the redundant coordinate as the triangilation in 3d fails for a
% % plane surface
% Mnp1=Mnp(:,[1,3]);
% DT = delaunayTriangulation(Mnp1);
% nodeCoord1=nodeCoord(:,[1,3]);

%% % plot tria mesh
% % this block of code did not work for the 3D case as point location is
% % notcompatible with 3d tri
% tri=load('tria.txt');
% tri=tri+ones(size(tri));
% coOrdinate= load('half.xyz');
% trimesh(tri,coOrdinate(:,1),coOrdinate(:,2),coOrdinate(:,3));
% axis equal;
%
% DT = triangulation(tri,coOrdinate);

%%
% [ti,Ba] = pointLocation(DT,nodeCoord) ; % Nan correspnds to the points outside the grid
% Ba(isnan(Ba))=0;  % Ba is the barycentric coordinate
% ti(isnan(ti))=0;
%
% eigVec=zeros(size(nodeCoord,1),size(V,2));
%
% for i=1:length(ti)
%
%     elementID=ti(i);
%     if elementID==0   % if node is not inside the triangulation then consider the nearest node
%         nearestNode=dsearchn(Mnp,nodeCoord(i,:));
%         eigVec(i,:)=V(nearestNode,:);
%     else              % interpolation of eigen values according to wikle and cressie
%         tria=DT.ConnectivityList(elementID,:);
%         baryCoord=Ba(i,:);
%         eigVec(i,:)=baryCoord(1).*V(tria(1),:)+baryCoord(2).*V(tria(2),:)+baryCoord(3).*V(tria(3),:);
%
%     end
%
% end

%% loading interpolated eigen vectors
eigVec=load('interpEigVecInner3_400.mat');%interEigVecs.mat');%
eigVec=eigVec.interpEigVec;
%% klaman recursion

% at time zero
a00=pc(1,:)';
p00=eye(size(V,2));
ytp1=zeros(size(nodeCoord,1),size(devCenterd,1));
ytp1=ytp1';
ytp=zeros(size(ytp1));
zt=devCenterd';

for i=1:5%size(zt,2)
	
	
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
	
	ytp(i,:)=att'*eigVec';   % taking transpose to maintain convention of deviation matrix
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


for i=1:2
	set(0, 'currentfigure', fig1)
	plotAxis1=subplot(2,2,i);
	contourDomainPlot(fem,domainID,devCenterd1(i,:),plotAxis1)
	ax=gca;
	ax.Clipping='off';
	%     caxis([-2 2]);
	
	set(0, 'currentfigure', fig2)
	plotAxis2=subplot(2,2,i);
	contourDomainPlot(fem,domainID,ytp(i,:),plotAxis2)
	ax=gca;
	ax.Clipping='off';
	hold all
	plot3(nodeCoord(iMnp,1),nodeCoord(iMnp,2),nodeCoord(iMnp,3),'*');
	%     caxis([-2 2]);
	
	set(0, 'currentfigure', fig3)
	plotAxis3=subplot(2,2,i);
	contourDomainPlot(fem,domainID,ytp1(i,:),plotAxis3)
	ax=gca;
	ax.Clipping='off';
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
stringTitle=sprintf('RMSE for prediction in time t after measuring %d points',length (iMnp));
title(stringTitle);
xlabel('instance');
ylabel('RMSE in mmm');

figure;
rmseT=sqrt(sum((bsxfun(@minus,devCenterd1(2:end,:),ytp1(1:end-1,:))).^2,2)./size(devCenterd1,2));
plot(rmseT);
stringTitle=sprintf('RMSE for prediction for time (t+1) after measuring %d points in time t',length (iMnp));
title(stringTitle);
xlabel('instance');
ylabel('RMSE in mmm');


