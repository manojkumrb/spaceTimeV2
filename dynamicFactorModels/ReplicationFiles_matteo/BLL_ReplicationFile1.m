%% BLL_ReplicationFile1 - Impulse Response Function for Monetary Policy Shock

clear all; close all; clc; 
ML_graph_options
sel=2; trans=2;
r=7; q=3; d=2; s=200; p=2; det=1;                                           % parameter set up
[Y, Label, Name, Dates, cd] = BLL_ReadData(3,1,sel);                        % Load data in levels
[Data,~,~,~, codeData] = BLL_ReadData(2,1,sel);                             % Load Data in first differences

%% Set Parameters for Structural Analysis %%
boot=1000;                                                                 	% number of bootstrap
idvar=[31 7 65];                                                            % variables on which restrictions are imposed
LRR=0; SRR=[1 2;1 3;2 3];                                                   % choleski

%% Structural Estimation %%
[C2,e] = BLL_SLRImp(Y,cd,q,d,r,p,s,idvar,SRR,LRR);                          % IRF - BLL (VECM)
[CC2, eflag2] = BLL_SLRBoot(Y,cd,q,d,r,p,s,idvar,SRR,LRR,boot);             % Bands - BLL (VECM)
[C3,sh2] = BLL_SLRImpVAR(Y,cd,q,r,p,s,idvar,SRR,LRR);                       % IRF - BLL (VAR)
[CC3, eflag3] = BLL_SLRBootVAR(Y,cd,q,r,p,s,idvar,SRR,LRR,boot);            % Bands - BLL (VAR)

%% Rearrange and Normalize Impulse Response Functions %%
qq=3;
IRF{1,1}=squeeze(C2(:,qq,:))';
IRF{1,2}=squeeze(C3(:,qq,:))';
for bb=1:boot;
    IRFboot{1,1}(:,:,bb)=squeeze(CC2(:,qq,:,bb))';
    IRFboot{1,2}(:,:,bb)=squeeze(CC3(:,qq,:,bb))';
end


for mm=1:2;    
    IRF{1,mm}=ML_IRFnormalize(IRF{1,mm},.5,67);                   % Monetary Policy Shock
    for bb=1:boot;        
        IRFboot{1,mm}(:,:,bb)=ML_IRFnormalize(IRFboot{1,mm}(:,:,bb),.5,67);
    end
end


%% Confidence bands %%
xi=68; up=.5*(100+xi); down=.5*(100-xi);
for jj=1:2;       
    first=mean(IRFboot{1,jj},3)';
    second=mean(IRFboot{1,jj}.^2,3)';
    stde=sqrt(second-first.^2)';
    IRFu{1,jj} = IRF{1,jj} + stde;
    IRFd{1,jj} = IRF{1,jj} - stde;
end


%% Graphs %%
ID=[31 7 65 34];                                                     % Variables of Interests
Variable={'Gross Domestic Product','Consumer Price Index','Federal Funds Rate','Residential Investments'};

close all;
s2=81;
for ii=1:4;
    a=min(0,ML_min([IRFd{1,2}(1:s2,ID(ii)) IRFd{1,1}(1:s2,ID(ii))]));
    b=max(0,ML_min([IRFu{1,2}(1:s2,ID(ii)) IRFu{1,1}(1:s2,ID(ii))],2));
    c=round(10*(b-a)/10)/10;
    axes('Parent',figure,'FontSize',10); ML_FigureSize(1); hold on;    
    ha=area(0:s2-1,[IRFd{1,2}(1:s2,ID(ii)) IRFu{1,2}(1:s2,ID(ii))-IRFd{1,2}(1:s2,ID(ii))],'linestyle','none');
    set(ha(1), 'FaceColor', 'none'); set(ha(2), 'FaceColor', [0.8 0.8 0.8])        
    plot(0:s2-1,zeros(s2,1),'k');
    plot(0:s2-1,IRF{1,2}(1:s2,ID(ii)),'color',[.5 .5 .5],'linewidth',3); 
    plot(0:s2-1,IRFu{1,1}(1:s2,ID(ii)),'k--',0:s2-1,IRFd{1,2}(1:s2,ID(ii)),'k--','linewidth',2);
    plot(0:s2-1,IRF{1,1}(1:s2,ID(ii)),'k-','linewidth',3); 
    hold off; axis tight; box on
    set(gca,'Xtick',0:8:80,'Xticklabel',0:2:20,'Layer','top')
    set(gca,'Ytick',ML_ytick(a,b,c))
    xlabel('Years After the shock','fontangle','italic','fontsize',12);
    title(Variable{ii},'fontangle','italic','fontsize',18);           
    axis([xlim+eps a-c/10 b+c/10])
end