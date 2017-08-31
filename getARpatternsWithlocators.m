% simulate variation patterns for a component using FEM. case: flat plate
% getting deviaiton pattens with slots holes and clamps; the deviations are
% auto correlated
close all;clear;

run 'C:\Users\babu_m\code\source\initFiles.m';
fem=femInit(source);

fileName{1}=strcat(inpFiles,'\top_hat.inp');%\hinge.inp');%\halo.inp');%
fem=importMesh(fem,fileName{1});%'top_hat.inp'); % has provisions for domain separation within code
totalNodes=size(fem.xMesh.Node.Coordinate,1);
% fem=importMultiMesh(fem, fileName);
fem=femPreProcessing(fem);
domainID=2;
idPart=domainID;
normal=fem.xMesh.Node.Normal;
coordi=fem.xMesh.Node.Coordinate;
nodeIdDomain=fem.Domain(idPart).Node;
nodeCoord=fem.xMesh.Node.Coordinate(nodeIdDomain,:); % coordinates for the domain
% ar parameters
a0=0;
a1=0.6;
sigma=1;
nPatterns=9;
% locators' details
nodeIdNc=nodeIdDomain([625,630]);
nodeIdClamps=nodeIdDomain([1283,1286,26,247]);
nodeIdSlots=nodeIdDomain(796);
nodeIdPin=nodeIdDomain(422);
dir1n=[0,1,0];
dofClamps=2;
dir1s=[1,0,0];
dir2s=[0,0,1];
dir1p=[1,0,0];
dir2p=[0,0,1];
%get ar deviations
devNc=getAR1Dev(a0,a1,sigma,length(nodeIdNc),nPatterns);
devClamps=getAR1Dev(a0,a1,sigma,length(nodeIdClamps),nPatterns);

for i=1:nPatterns
    
    fem=setNCblocks(fem,nodeIdNc,dir1n,domainID,devNc(i,:),[]);
    fem=setClamps(fem,nodeIdClamps,dofClamps,devClamps(i,:));
    fem=setSlots(fem,nodeIdSlots,dir1s,dir2s,domainID,[]);
    fem=setPin(fem,nodeIdPin,dir1p,dir2p,domainID,[]);
    
    fem=femReset(fem);
    fem=femRefresh(fem);
    fem=femSolve(fem);
    
    % node deviations are in fem.Sol.U, u,v,w,rot1,rot2,rrot3 sequentially for all nodes
    femDev=fem.Sol.U(2:6:end);% 1 for u, 2 for v,3 for w etc..
    %%%the deviaiton of interest%%%%%%%%%%%%%%%%%%%%%%
    femDevDomain(i,:)=femDev(nodeIdDomain);
    
    %% Plotting the boundary constraints
    figure
    ax=gca;
    hold all
    fem.Post.Options.ParentAxes=ax;
    fem.Post.Contour.Domain=domainID;
    fem.Post.Contour.ContourVariable='v';
    fem.Post.Contour.Resolution=1;
    fem.Post.Options.SymbolSize=3;
    fem.Post.Contour.MaxRange=5;
    fem.Post.Contour.MinRange=-5;
    contourPlot(fem);
    bilateralNodeBcPlot(fem)
    unilateralBcPlot(fem)
    pinslotBcPlot(fem)
    caxis([-2 2]);
        fem.Post.Options.LengthAxis=0;
    %     loadPlot(fem)
    axis off
end



%% plotting time series of deviations of clamps and locators
figure
plot(devNc);hold all; plot(devClamps);
xlabel('n^{th} time instance');
ylabel('deviation in mm');
legend('NC Block 1','NC Block 2','clamp 1','clamp 2','clamp 3','clamp 4');

% femDevArFlatPlate=femDevDomain;


