% comparing deterministinc and adaptive measurement
% reads the result from all possible determisnistic sequences and then
% finds the sequence with minimum rmse to compare with the adaptive
% measurements
% variable delta stores the best sequence from all possible cases

close all;clear all;
run 'C:\Users\babu_m\Documents\GitHub\source\initFiles.m';
path(path,'\\megara\dlmr\Manoj\all data sensitivity inner'); % all data files in megara

load('innerEigWithkrigYtp12batch35basis_1.0W1_0.0W2.mat'); % contains structure seqPred %


%% Rmse diff AFTER EACH MEASUREMET AVERAGED FOR ALL PATTERNS

numberCompleteMeasure   =10;
maxSnap                 =9;
nPattern                =length(seqPred);
totalInstances          = nPattern-(numberCompleteMeasure);

regionsMes(maxSnap).Seq(1).Pattern(1).Dev=[];

%% rmse for a given sequence(j) averaged over all parts(variation pattern(i)), compared
%with average of rmse of adaptive measurements after same no of
%measurements for all parts(variation pattern)
for nMes=1:maxSnap
    fileName=sprintf('detAllSeq%dReg.mat',nMes);
    tempvar=load(fileName);
    regionsMes(nMes).Seq=tempvar.regionsMes.Seq;
    
    delta(nMes).Min=inf;
    delta(nMes).MinSeq=[];
    
    delta(nMes).Max=0;
    delta(nMes).MaxSeq=[];
    
    for j=1:length(regionsMes(nMes).Seq)
        rmseD=0;
        rmseA=0;
        
        allMesReg=regionsMes(nMes).Seq(j).Regions;
        nMesReg=length(allMesReg);
        for i=numberCompleteMeasure+1:totalInstances % i for each part
            
            % adding Rmse
            rmseD=rmseD+regionsMes(nMes).Seq(j).Pattern(i).Rmse;
            rmseA=rmseA+seqPred(i).Snap(nMesReg).RmseT(nMesReg);
        end
        delta(nMes).Seq(j)=((rmseD-rmseA)/rmseD)*100;
        
        if delta(nMes).Seq(j)<delta(nMes).Min
            delta(nMes).Min=delta(nMes).Seq(j);
            delta(nMes).RmseDetMin=rmseD;
            delta(nMes).RmseA=rmseA;
            delta(nMes).MinSeq=allMesReg;
            
        end
        if delta(nMes).Seq(j)>delta(nMes).Max
            delta(nMes).Max=delta(nMes).Seq(j);
            delta(nMes).RmseDetMax=rmseD;
            delta(nMes).RmseA=rmseA;
            delta(nMes).MaxSeq=allMesReg;
        end
    end
    
end

%% plotting 
for i=1:maxSnap
    rmseAvgDetMin(i)=delta(i).RmseDetMin;
    rmseAvgDetMax(i)=delta(i).RmseDetMax;
    rmseAvgAdap(i)=delta(i).RmseA;
end

rmseAvgDetMin	= rmseAvgDetMin-rmseAvgDetMin(maxSnap).*ones(size(rmseAvgDetMin));
rmseAvgDetMax	= rmseAvgDetMax-rmseAvgDetMax(maxSnap).*ones(size(rmseAvgDetMax));
rmseAvgAdap		= rmseAvgAdap-rmseAvgAdap(maxSnap).*ones(size(rmseAvgAdap));

% standard error claculation
rmseSE=zeros(totalInstances,maxSnap);
count=1;
for i=numberCompleteMeasure+1:nPattern    
    rmseSE(count,:)=seqPred(i).Snap(maxSnap).RmseT;
    count=count+1;
end
stdError=std(rmseSE)./sqrt(size(rmseSE,1));
stdError=stdError.*1.96;
% average
rmseAvgDetMin=rmseAvgDetMin./totalInstances;
rmseAvgAdap=rmseAvgAdap./totalInstances;
rmseAvgDetMax=rmseAvgDetMax./totalInstances;

% plotting standard error
upperStdErr=rmseAvgAdap+stdError;
lowerStdErr=rmseAvgAdap-stdError;
fillRmseYAdap=[upperStdErr,fliplr(lowerStdErr)]';

fig=figure;

fillX=[1:length(rmseAvgDetMax),fliplr(1:length(rmseAvgDetMin))]';
fillYDet=[rmseAvgDetMax,fliplr(rmseAvgDetMin)]';

patchHandle=patch(fillX,fillYDet,[224,243,219]./256);
patchHandle.FaceAlpha=0.4;
patchHandle.EdgeColor=[255,255,255]./256; % now white
hold all
patchHandle=patch(fillX,fillRmseYAdap,[254,178,76]./256);
patchHandle.FaceAlpha=0.4;
patchHandle.EdgeColor=[255,255,255]./256;

% Plotting means, max and min
hold all
hh=plot(1:length(rmseAvgDetMin),rmseAvgDetMin,...
    1:length(rmseAvgDetMax),rmseAvgDetMax,...
    1:length(rmseAvgAdap),rmseAvgAdap);

fig=changeAxesLooks(fig,{'Average RMSE of Best deterministic prediciton','compared to adaptive prediction'},...
    'Number of measurements','RMSE in mm');

set(hh,'linewidth',1.5);

legendflex(hh, {'Best deterministic prediciton','Worst deterministic prediciton','Adaptive prediction'}, ...
    'anchor', {'ne','ne'}, ...
    'buffer', [0 0], ...
    'nrow',   3,...
    'fontsize',12, ...
    'xscale', 0.66, ...
    'box','on');


print('-r400','detAdpCompHinge','-dpng');


