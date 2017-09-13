% simulate variation patterns for a component using FEM. case: hinge
% getting deviaiton pattens with slots holes and clamps; the deviations are
% auto correlated
close all;clear;

run 'C:\Users\babu_m\Documents\GitHub\source\initFiles.m';
fem=femInit();

fileName{1}=strcat(inpFiles,'\Remesh_hinge.inp');%\top_hat.inp');%\hinge.inp');%\Remesh_hinge.inp');%\halo.inp');%
fem=importMesh(fem,fileName{1});%'top_hat.inp'); %use only for top hat
% fem=importMultiMesh(fem, fileName);%has provisions for domain separation 
% within code
fem=femPreProcessing(fem);
domainID=2;
idPart=domainID;
normal=fem.xMesh.Node.Normal;
coordi=fem.xMesh.Node.Coordinate;
nodeIdDomain=fem.Domain(idPart).Node;
nodeCoord=fem.xMesh.Node.Coordinate(nodeIdDomain,:); % coordinates for the domain

% %% finding mesh nodeID

% meshplotAxisDefined(fem,domainID); figure=gcf; axiss=gca;
% [nodeIDout,nodeIDoutRect]=guiSelection(figure,axiss,nodeCoord);

% ar parameters
a0=0;
a1=0.6;
sigma=1;
nPatterns=5;

%% locators' details (getting global node coordinates)
% % For full mesh
% %>>>>>>>>the meshimport function is the one used only for % top hat<<<<<<<<
% nodeIdNc=nodeIdDomain([10841,6119,661,10749,1392,8063]);
% nodeIdClamps=nodeIdDomain([12026,10396,11394,3841]);
% nodeIdSlots=nodeIdDomain(2232);
% nodeIdPin=nodeIdDomain(5008);

% For coarse mesh
nodeIdNc=nodeIdDomain([1108,686,675,416,387,238]);
nodeIdClamps=nodeIdDomain([1128,1124,258,255]);
nodeIdSlots=nodeIdDomain(817);
nodeIdPin=nodeIdDomain(376);

dir1n=normal(nodeIdNc,:);%[0,1,0];
dir1n(1:2:end,:)=dir1n(1:2:end,:).*-1;%[0,1,0]; % reversing alternate clamp direction
dirC=normal(nodeIdClamps,:);
slotNormal=normal(nodeIdSlots,:);
pinNormal=normal(nodeIdPin,:);

dir1s=zeros(length(nodeIdSlots),3);
dir2s=zeros(length(nodeIdSlots),3);
dir1p=zeros(length(nodeIdPin),3);
dir2p=zeros(length(nodeIdPin),3);

for i=1:length(nodeIdSlots)
    basis=null(slotNormal(i,:));
    dir1s(i,:)=basis(:,1);
    dir2s(i,:)=basis(:,2);
end

for i=1:length(nodeIdPin)
    basis=null(pinNormal(i,:));
    dir1p(i,:)=basis(:,1);
    dir2p(i,:)=basis(:,2);
end

%get ar deviations
devNc=getAR1Dev(a0,a1,sigma,length(nodeIdNc),nPatterns);
devClamps=getAR1Dev(a0,a1,sigma,length(nodeIdClamps),nPatterns);
femDevDomain=zeros(nPatterns,length(nodeIdDomain));

for i=1:nPatterns
    
    fem=setNCblocks(fem,nodeIdNc,dir1n,domainID,devNc(i,:),[]);
    fem=setClampsNormalTra(fem,nodeIdClamps,dirC,devClamps(i,:));
    %     fem=setClamps(fem,nodeIdClamps,2,devClamps(i,:));
    fem=setSlots(fem,nodeIdSlots,dir1s,dir2s,domainID,[]);
    fem=setPin(fem,nodeIdPin,dir1p,dir2p,domainID,[]);
    fem.Options.Solver.MaxIter=15;
    fem.Options.Solver.Eps=1e-3;
    
    fem=femReset(fem);
    fem=femRefresh(fem);
    fem=femSolve(fem);
    
    % node deviations are in fem.Sol.U, u,v,w,rot1,rot2,rrot3 sequentially for all nodes
    femDev(:,1)=fem.Sol.U(1:6:end);% 1 for u, 2 for v,3 for w etc..
    femDev(:,2)=fem.Sol.U(2:6:end);% 1 for u, 2 for v,3 for w etc..
    femDev(:,3)=fem.Sol.U(3:6:end);% 1 for u, 2 for v,3 for w etc..
    %%%the deviaiton of interest%%%%%%%%%%%%%%%%%%%%%%
    femDevDomain(i,:)=sqrt(sum(femDev(nodeIdDomain,:).^2,2));
%     femDevDomain(i,:)=femDevDomain(i,nodeIdDomain);
    %% Plotting the boundary constraints
%         figure;
%         ax=gca;
%         hold all
%         fem.Post.Options.ParentAxes=ax;
%         fem.Post.Contour.Domain=domainID;
%         fem.Post.Contour.ContourVariable='user';
%         fem.Sol.UserExp=femDevDomain(i,:);
%         fem.Post.Contour.Resolution=1;
%         fem.Post.Options.SymbolSize=5;
%         fem.Post.Contour.MaxRange=5;
%         fem.Post.Contour.MinRange=-5;
%         contourPlot(fem);
%         bilateralNodeBcPlot(fem)
%         unilateralBcPlot(fem)
%         pinslotBcPlot(fem)
%         caxis([-2 2]);
%         fem.Post.Options.LengthAxis=0;
%         loadPlot(fem)
%         axis off
contourDomainPlot(fem,1,femDevDomain(i,:),0)
contourcmap('parula',0.25);
end



%% plotting time series of deviations of clamps and locators

myTtile='Auto correlated locators'' variation'; % two consequtive apostrophe for printing one apostrophe
myXlabel='Time (instance)';
myYlabel='Deviation in mm';


h=figure('Units','centimeters');%210 × 297 millimeters 
% h.Position=[5,10,15,10];
h.PaperPositionMode='auto';

[ha, pos] = tight_subplot(2, 1, 0.02, [0.2 0.1], [0.1 0.05]);
markers = {'o','s','d','^','v','x','+','*','.','>','<','p','h'};

ncAxis=ha(1);
ncLine=plot(ncAxis,devNc);


clampAxis=ha(2);
clampLine=plot(clampAxis,devClamps);



cmapNC = cbrewer('seq','Reds',2*length(ncLine));
cmapNC=flipud(cmapNC);

% cmapNC = cbrewer('seq','Greys',length(clampLine));
cmapClamp = cbrewer('seq','Blues',2*length(clampLine));
cmapClamp=flipud(cmapClamp);

legNcCell=cell(length(ncLine),1);
for i=1:length(ncLine)
   ncLine(i).LineWidth=1.5;
   ncLine(i).Color=cmapNC(i,:); 
   ncLine(i).Marker=markers{i};
   ncLine(i).MarkerFaceColor=cmapNC(i,:); 
   ncLine(i).MarkerEdgeColor='none';
   legNcCell{i}=sprintf('%i',i);
end

legClampCell=cell(length(clampLine),1);
for i=1:length(clampLine)
   clampLine(i).LineWidth=1.5;
   clampLine(i).Marker=markers{i};
   clampLine(i).Color=cmapClamp(i,:); 
   clampLine(i).MarkerFaceColor=cmapClamp(i,:);
   legClampCell{i}=sprintf('%i',i);

end

axes(clampAxis)                                             % Makes clamp axis current
myAxis=h.CurrentAxes;                                       % queries the axis from the figure
set(myAxis, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.01 .01] , ...
  'xTick'       , 1:1000 , ...
  'xTicklabel'  , ['0';'1';'2';'3';'4';'5';'6'] , ...
  'XMinorTick'  , 'off'     , ...
  'YMinorTick'  , 'off'     , ...
  'YGrid'       , 'on'      , ...
  'XGrid'       , 'off'      , ...
  'Ylim'        , [-4 4]     , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'YTick'       , -4:1:4, ...
  'LineWidth'   , 1,...
  'FontName'    , 'Times New Roman',...
  'fontsize'    , 12        );

axes(ncAxis)                                                % Makes NCclamp axis current
myAxis=h.CurrentAxes;                                       % queries the axis from the figure
set(myAxis, ...
  'Box'         , 'off'     , ...
  'xTick'       , []        , ...
  'YTick'       , -4:1:4, ...
  'TickDir'     , 'out'     , ...
  'Ylim'        , [-4 4]     , ...
  'TickLength'  , [.01 .01] , ...
  'XMinorTick'  , 'off'     , ...
  'YMinorTick'  , 'off'     , ...
  'XColor'      , [1 1 1 ], ...
  'YColor'      , [.3 .3 .3], ...
  'YGrid'       , 'on'      , ...
  'XGrid'       , 'on'      , ...
  'LineWidth'   , 1,...
  'FontName'    , 'Times New Roman',...
  'fontsize'    , 12        );

[hl(2).leg, hl(2).obj, hl(2).hout, hl(2).mout] = ...
legendflex(ncLine, legNcCell, ...
'anchor', {'ne','se'}, ...
'buffer', [0 -20], ...
'nrow',   1,...
'fontsize',10, ...
'xscale', 0.66, ...
'box','off',...
'title', 'NC Block');

[hl(3).leg, hl(3).obj, hl(3).hout, hl(3).mout] = ...
legendflex(clampLine, legClampCell, ...% 'ref', hl(2).leg, ...
'anchor'    ,{'ne','se'}, ...
'buffer'    ,[0 0]      , ...
'nrow'      ,1          ,...
'fontsize'  ,10         , ...
'box'       ,'off'       ,...
'title'     ,'Clamp'   );

figtitle(myTtile,...
    'FontSize'   , 16          , ...
    'FontWeight' , 'bold'     );

[~,hXLabel]=suplabel(myXlabel,'x',[.1 .18 .84 .84]);
[~,hYLabel]=suplabel(myYlabel,'y',[.1 .1 .84 .84]);

set([hXLabel, hYLabel]  , ...
'FontSize'   , 16       , ...
'FontWeight' , 'bold'   );

% print(h,myTtile,'-dpng')


