% BLL_ReplicationFile2 - Impulse Responses for Supply shock 

clear all; close all; clc; 
ML_graph_options
sel=2; trans=2;
r=7; q=3; d=2; s=200; p=2; det=1;                                           % parameter set up
[Y, Label, Name, Dates, cd] = BLL_ReadData(3,1,sel);                        % Load data in levels
[Data,~,~,~, codeData] = BLL_ReadData(2,1,sel);                             % Load Data in first differences
boot=1000;                                                                  % number of bootstrap
idvar=[31 33 44];                                                           % variables on which restrictions are imposed
LRR=[1 2;1 3]; SRR=[1 3];                                                   % Short and Long Run Restrictions
Stima=0;

%%% --------------------- %%%
%%% Structural Estimation %%%
%%% --------------------- %%%
[C2,e] = BLL_SLRImp(Y,cd,q,d,r,p,s,idvar,SRR,LRR);                          % IRF - BLL (VECM)
CC2 = BLL_SLRBoot(Y,cd,q,d,r,p,s,idvar,SRR,LRR,boot);                       % Bands - BLL (VECM)
[C3,sh2] = BLL_SLRImpVAR(Y,cd,q,r,p,s,idvar,SRR,LRR);                       % IRF - BLL (VAR)
CC3 = BLL_SLRBootVAR(Y,cd,q,r,p,s,idvar,SRR,LRR,boot);                      % Bands - BLL (VAR)

%%% -------------------------------------------------- %%%
%%% Rearrange and Normalize Impulse Response Functions %%%
%%% -------------------------------------------------- %%%
IRF{1,1}=squeeze(C2(:,1,:))';
IRF{1,2}=squeeze(C3(:,1,:))';
for bb=1:boot;    
    IRFboot{1,1}(:,:,bb)=squeeze(CC2(:,1,:,bb))';
    IRFboot{1,2}(:,:,bb)=squeeze(CC3(:,1,:,bb))';
end


for mm=1:2;
    IRF{1,mm}=ML_IRFnormalize(IRF{1,mm},.25,31);    
    for bb=1:boot;
        IRFboot{1,mm}(:,:,bb)=ML_IRFnormalize(IRFboot{1,mm}(:,:,bb),.25,31);
    end
end

%%% ---------------- %%%
%%% Confidence bands %%%
%%% ---------------- %%%
xi=68; up=.5*(100+xi); down=.5*(100-xi);
for jj=1:2;
    temp=IRFboot{1,jj}-repmat(IRF{1,jj},[1 1 boot]);
    IRFu{1,jj}=IRF{1,jj}-prctile(temp,down,3); IRFd{1,jj}=IRF{1,jj}-prctile(temp,up,3);
end


%%% -------------------------- %%%
%%% Plotting Impulse Responses %%%
%%% -------------------------- %%%
ID = [31 34 35 40:42];
Variabile={'Gross Domestic Product';'Residential Investment'; 'Nonresidential Investment';
    'Consumption: Nondurable Goods';'Consumption: Services';'Consumption: Durable Goods'};
for ii=1:6;      
    a=min(0,ML_min([IRFd{1,2}(1:s,ID(ii)) IRFd{1,1}(1:s,ID(ii))]));
    b=max(0,ML_min([IRFu{1,2}(1:s,ID(ii)) IRFu{1,1}(1:s,ID(ii))],2));
    c=round(10*(b-a)/10)/10;
    axes('Parent',figure,'FontSize',10); ML_FigureSize(1); hold on;
    plot(0:s-1,zeros(s,1),'k');
    ha=area(0:s-1,[IRFd{1,2}(1:s,ID(ii)) IRFu{1,2}(1:s,ID(ii))-IRFd{1,2}(1:s,ID(ii))],'linestyle','none');
    set(ha(1), 'FaceColor', 'none'); set(ha(2), 'FaceColor', [0.8 0.8 0.8])    
    plot(0:s-1,IRF{1,2}(1:s,ID(ii)),'Color',[0.5 0.5 0.5],'linewidth',3); 
    plot(0:s-1,IRFu{1,1}(1:s,ID(ii)),'k--',0:s-1,IRFd{1,1}(1:s,ID(ii)),'k--','linewidth',2);
    plot(0:s-1,IRF{1,1}(1:s,ID(ii)),'k-','linewidth',3);          
    axis tight; hold off; box on;
    set(gca,'Xtick',0:20:200,'Xticklabel',0:5:50,'Layer','top')
    xlabel('Years After the Shock','fontangle','italic','fontsize',12);
    title(Variabile{ii},'fontangle','italic','fontsize',18);        
    set(gca,'Ytick',ML_ytick(a,b,c))
    axis([xlim+eps a-c/10 b+c/10])     
end

