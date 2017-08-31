% comparing deterministinc and adaptive measurement


close all;clear all;
run 'C:\Users\babu_m\code\source\initFiles.m';

load('hingeDevArSimLocCoarse.mat'); % contains variable femDevDomain
load('seqPedStructHingeCoarseWithRmseTp1_25Eig.mat'); % contains structure seqPred %
load('detSeqHingeCoarse1to5_40rep.mat'); % contains structure seqPred %

zt=femDevDomain;


%% deterministic pattern search
numberCompleteMeasure=10;
regionMes(4).Seq(15).Pattern(1).Dev=[];

for nMes=1:5
    c = nchoosek(1:6,nMes);
    delta(nMes).Min=inf;
    delta(nMes).MinSeq=[];
    
    for j=1:size(c,1)
        

        rmseD=0;
        rmseA=0;
        allMesReg=regionMes(nMes).Seq(j).Regions;
        nMesReg=length(allMesReg);
        for i=numberCompleteMeasure+1:size(zt,1) % i for each part
            
            % saving predictions
            ytt(i,:)=regionMes(nMes).Seq(j).Pattern(i).Dev;
            VarYtt(i,:)=regionMes(nMes).Seq(j).Pattern(i).Var;
            rmseD=rmseD+regionMes(nMes).Seq(j).Pattern(i).Rmse;
            rmseA=rmseA+seqPred(i).Snap(nMesReg).RmseT(nMesReg);
        end
        delta(nMes).Seq(j)=((rmseD-rmseA)/rmseD)*100;
        
        if delta(nMes).Seq(j)<delta(nMes).Min
            delta(nMes).Min=delta(nMes).Seq(j);
            delta(nMes).MinSeq=allMesReg;
        end
    end
    
    bestSeqComp(nMes)=delta(nMes).Min;
end

fig=figure;
hh=plot(bestSeqComp(2:end));
fig=changeAxesLooks(fig,{'Average improvement in RMSE with adaptive measurement','compared to best fixed measurement sequence'},...
    'Number of measurements','Percentage improvement');
hh.LineWidth=1.5;


