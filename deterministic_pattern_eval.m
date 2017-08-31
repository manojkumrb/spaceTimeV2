% get performance of combination of all possible deterministic sequences and saving the values for
% further evaluation.
close all;
clear ;
% dbstop if error;
run 'C:\Users\babu_m\Documents\GitHub\source\initFiles.m';
path(path,'\\megara\dlmr\Manoj\all data sensitivity inner'); % all data files in megara

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
load('optimisedHypParmHinge.mat'); %contains struct hyp in gp function format(i.e. hyp.cov=log(parameters)) for matern covariance
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

%% deterministic pattern search

regionMes(4).Seq(15).Pattern(1).Dev=[];

for nMes=1:5
	
	% creating all possible permutations
	c = nchoosek(1:6,nMes);
	count=1;
	cTotal=zeros(size(c,1)*factorial(nMes),nMes);
	
	for j=1:size(c,1)
			cc=perms(c(j,:));
		
		for jj=1:size(cc,1)
			
			cTotal(count,:)=cc(jj,:);
			count=count+1;
			
		end
	end
	
	for kk=1:size(cTotal,1)
		
		for i=numberCompleteMeasure+1:size(zt,1) % i for each part
			
			allMesReg=cTotal(kk,:); % all the measured regions
			[tempV,tempR,tempEigVec,tempZ]=getRegionSysMat(nodeIDoutRect,allMesReg,V,R,eigVec,zt(i,:));
			% updating the state variables after partial measurement
			[tempKt,att,ptt]=getKalmanStatett(pttm1,attm1,tempEigVec,tempR,tempV,tempZ');
			% predicting deviations for the whole part using updated state variables
			[ytt(i,:),VarYtt(i,:),ytp1(i,:),VarYtp1(i,:)]=getPredictions(eigVec,att,ptt,H,V,R,covErr);
			
			% RMSE calculation for the step
			diff=bsxfun(@minus,devCenterd(i,:),ytt(i,:));
			rmseT=sqrt(sum((diff).^2)/size(devCenterd,2));  % rmse without noisy input
			
			% saving predictions
			regionMes(nMes).Seq(kk).Regions=allMesReg;
			regionMes(nMes).Seq(kk).Pattern(i).Dev=ytt(i,:);
			regionMes(nMes).Seq(kk).Pattern(i).Var=VarYtt(i,:);
			regionMes(nMes).Seq(kk).Pattern(i).Rmse=rmseT;
		end
		fprintf('%d of %d sequences completed with %d regions \n',kk,size(cTotal,1),nMes);
		
		
	end
	
	
	
end
