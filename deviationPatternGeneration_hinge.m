% Old version not used 

%simulate variation patterns for a component using FEM. case: flat plate
% 12 nodes selected and randomly alotted deviations of -2 or 0 or 2 with equal probability
close all;clear all;

run 'C:\Users\babu_m\code\source\initFiles.m';
fem=femInit(source);

fileName{1}=strcat(inpFiles,'\hinge.inp');%\top_hat.inp');%\halo.inp');%
% dev=-ones(size(fem.xMesh.Node.Coordinate,1),1)*1.8;
% fem=createVariationalMesh(fem, dev); fem=femPreProcessing(fem);
fem=importMesh(fem,fileName{1});%'top_hat.inp'); % has provisions for domain separation within code
totalNodes=size(fem.xMesh.Node.Coordinate,1);
% fem=importMultiMesh(fem, fileName);
fem=femPreProcessing(fem);
domainID=1;
idPart=domainID;
normal=fem.xMesh.Node.Normal;
coordi=fem.xMesh.Node.Coordinate;
nodeIdDomain=fem.Domain(idPart).Node;
nodeCoord=fem.xMesh.Node.Coordinate(nodeIdDomain,:); % coordinates for the domain
nodeCoord2=nodeCoord(:,[1,3]); % make node 2d


%% finding mesh nodeID

meshplotAxisDefined(fem,domainID); figure=gcf; axiss=gca;
[nodeIDout,nodeIDoutRect]=guiSelection(figure,axiss,nodeCoord);

%% fem load addition

nPatterns=5;
mNodeId=nodeIdDomain([1283;1278;1286;797;796;810;421;422;437;26;32;247]);
mNodes= length(mNodeId);
devSet=[-2;0;2];
probSel=1/length(devSet); % setting uniform probability of selection

dev=zeros(nPatterns,mNodes);
femDevDomain=zeros(nPatterns,length(nodeIdDomain));
DevType='AR';%'random';%






for i=1:nPatterns    
    
    if strcmpi(DevType,'AR')
        a0=0;
        a1=0.4;
        sigma=1;
        
        if i==1            
            dev(i,:)=a0;            
        else
            
            dev(i,:)=a0.*ones(1,mNodes)+a1.*dev(i-1,:)+randn(1,mNodes).*sigma;
            setZeroDev=randi(mNodes,6,1);
            dev(i,setZeroDev)=0;
            
        end
        
    elseif strcmpi(DevType,'Random')
        
        for j=1:mNodes
            temp=rand();
            if temp <= probSel
                dev(i,j)=devSet(1);
                continue;
            elseif temp <= 2*probSel
                dev(i,j)=devSet(2);
                continue;
            else
                dev(i,j)=devSet(3);
            end
        end
    end
    
    fixedNodes=mNodeId(dev(i,:)==0);
    loadedNodes=mNodeId(dev(i,:)~=0);
    fem=setFixedNodes(fem, fixedNodes);
    fem=setNodeLoads(fem,loadedNodes,2*ones(mNodes),dev(i,:)); %[3185],[2],[0.5]);%
    fem=femRefresh(fem);
    fem=femSolve(fem);
    % node deviations are in fem.Sol.U, u,v,w,rot1,rot2,rrot3 sequentially for all nodes
    femDev=fem.Sol.U(2:6:end);% 1 for u, 2 for v,3 for w etc..
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%the deviaiton of interest%%%%%%%%%%%%%%%%%%%%%%
    femDevDomain(i,:)=femDev(nodeIdDomain);
    
    figure;
    ax=gca;
    hold all
    fem.Post.Options.ParentAxes=ax;
    fem.Post.Contour.Domain=domainID;
    fem.Post.Contour.ContourVariable='v';
    fem.Post.Contour.Resolution=1;
    fem.Post.Options.SymbolSize=10;
    fem.Post.Contour.MaxRange=inf;
    fem.Post.Contour.MinRange=-inf;
    contourPlot(fem);
    bilateralNodeBcPlot(fem)
    fem.Post.Options.LengthAxis=0.5;
    loadPlot(fem)
    axis off
    
end


% %% to check if femDevDomain has correct deviations
%     figure
%     ax=gca;
%     hold all
%     fem.Post.Options.ParentAxes=ax;
%     fem.Post.Contour.Domain=domainID;
%     fem.Post.Contour.ContourVariable='user';
%     fem.Sol.UserExp=femDev;
%     fem.Post.Contour.Resolution=1;
%     fem.Post.Options.SymbolSize=10;
%     fem.Post.Contour.MaxRange=inf;
%     fem.Post.Contour.MinRange=-inf;
%     contourPlot(fem);
%     bilateralNodeBcPlot(fem)
%     fem.Post.Options.LengthAxis=0.5;
%     loadPlot(fem)
%     axis off

