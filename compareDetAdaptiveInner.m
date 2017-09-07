% comparing deterministinc and adaptive measurement
% reads the result from all possible determisnistic sequences and then
% finds the sequence with minimum rmse to compare with the adaptive
% measurements
% variable delta stores the best sequence from all possible cases

close all;clear all;
run 'C:\Users\babu_m\Documents\GitHub\source\initFiles.m';
path(path,'\\megara\dlmr\Manoj\all data sensitivity inner'); % all data files in megara

load('innerEigWithkrigYtp12batch45basis_1.0W1_0.0W2.mat'); % contains structure seqPred %



% %% deterministic pattern search
% numberCompleteMeasure=10;
% regionsMes(9).Seq(1).Pattern(1).Dev=[];
% 
% for nMes=1:9
% 
% 	delta(nMes).Min=inf;
% 	delta(nMes).MinSeq=[];
% 	
% 	
% 	for j=1:length(regionsMes(nMes).Seq)
% 		
% 		
% 		rmseD=0;
% 		rmseA=0;
% 		allMesReg=regionsMes(nMes).Seq(j).Regions;
% 		nMesReg=length(allMesReg);
% 		
% 		for i=numberCompleteMeasure+1:50 % i for each part
% 			
% 			% saving predictions
% 			ytt(i,:)=regionsMes(nMes).Seq(j).Pattern(i).Dev;
% 			VarYtt(i,:)=regionsMes(nMes).Seq(j).Pattern(i).Var;
% 			rmseD=rmseD+regionsMes(nMes).Seq(j).Pattern(i).Rmse;
% 			rmseA=rmseA+seqPred(i).Snap(nMesReg).RmseT(nMesReg);
% 		end
% 		delta(nMes).Seq(j)=((rmseD-rmseA)/rmseD)*100;
% 		
% 		if delta(nMes).Seq(j)<delta(nMes).Min
% 			delta(nMes).Min=delta(nMes).Seq(j);
% 			delta(nMes).MinSeq=allMesReg;
% 		end
% 	end
% 	
% 	bestSeqComp(nMes)=delta(nMes).Min;
% end
% 
% fig=figure;
% hh=plot(bestSeqComp(1:end));
% fig=changeAxesLooks(fig,{'Average improvement in RMSE with adaptive measurement','compared to best fixed measurement sequence'},...
% 	'Number of measurements','Percentage improvement');
% hh.LineWidth=1.5;

%% Rmse diff AFTER EACH MEASUREMET AVERAGED FOR ALL PATTERNS

numberCompleteMeasure=10;
regionsMes(9).Seq(1).Pattern(1).Dev=[];


%% rmse for a given sequence(j) averaged over all parts(variation pattern(i)), compared
%with average of rmse of adaptive measurements after same no of
%measurements for all parts(variation pattern)
for nMes=1:9
	fileName=sprintf('detAllSeq%dReg.mat',nMes);
	tempvar=load(fileName);
	regionsMes(nMes).Seq=tempvar.regionsMes.Seq;
	delta(nMes).Min=inf;
	delta(nMes).MinSeq=[];
	
	for j=1:length(regionsMes(nMes).Seq)
		rmseD=0;
		rmseA=0;
		
		allMesReg=regionsMes(nMes).Seq(j).Regions;
		nMesReg=length(allMesReg);
		for i=numberCompleteMeasure+1:50 % i for each part
			
			% saving predictions
% 			ytt(i,:)=regionsMes(nMes).Seq(j).Pattern(i).Dev;
% 			VarYtt(i,:)=regionsMes(nMes).Seq(j).Pattern(i).Var;
			rmseD=rmseD+regionsMes(nMes).Seq(j).Pattern(i).Rmse;
			rmseA=rmseA+seqPred(i).Snap(nMesReg).RmseT(nMesReg);
		end
		delta(nMes).Seq(j)=((rmseD-rmseA)/rmseD)*100;
		
		if delta(nMes).Seq(j)<delta(nMes).Min
			delta(nMes).Min=delta(nMes).Seq(j);
			delta(nMes).RmseD=rmseD;
			delta(nMes).RmseA=rmseA;
			delta(nMes).MinSeq=allMesReg;
		end
	end
	
	bestSeqComp(nMes)=delta(nMes).Min;
end

for i=1:nMes
	rmseAvgDet(i)=delta(i).RmseD;
	rmseAvgAdap(i)=delta(i).RmseA;
end

nPattern=50;
totalInstances= nPattern-(numberCompleteMeasure);

% Plotting average RMSE
fig=figure;
hh=plot(1:length(rmseAvgDet),rmseAvgDet./totalInstances,1:length(rmseAvgAdap),rmseAvgAdap./totalInstances);
fig=changeAxesLooks(fig,{'Average RMSE of Best deterministic prediciton','compared to adaptive prediction'},...
	'Number of measurements','RMSE in mm');

set(hh,'linewidth',1.5);

legendflex(hh, {'Best deterministic prediciton','Adaptive prediction'}, ...
	'anchor', {'ne','ne'}, ...
	'buffer', [0 0], ...
	'nrow',   2,...
	'fontsize',12, ...
	'xscale', 0.66, ...
	'box','on');
print('-r400','detAdpCompHinge','-dpng');

% Plotting average RMSE improvement
fig=figure;
hh=plot(bestSeqComp(2:end));
fig=changeAxesLooks(fig,{'Average improvement in RMSE with adaptive measurement','compared to best fixed measurement sequence'},...
	'Number of measurements','Percentage improvement');
hh.LineWidth=1.5;


%% plotting rmse sensitivity separately
% 
% % nBasis		=40;%[5:5:60];
% numberCompleteMeasure=10;
% maxSnap=9;
% 
% % plot parameters
% markers = {'o','s','d','^','v','x','+','*','.','>','<','p','h','o','s','d','^','v','x'};
% cmapRmsePlot = cbrewer('qual','Set1',length(nBasis));
% 
% % reading data from file
% 
% for j=1:length(nTrainBatch)
% 	
% 	for w=1:length(infoWeight)
% 		
% 		for ww=1:length(distweight)
% 			
% 			for k=1:length(nBasis)
% 				
% 				fileString=sprintf('innerEigWithkrigYtp1%ibatch%ibasis_%.1fW1_%.1fW2.mat',nTrainBatch(j),nBasis(k),infoWeight(w),distweight(ww));
% 				seqPred=load (fileString);
% 				seqPred=seqPred.seqPred;
% 				
% 				for i=numberCompleteMeasure+1:length(seqPred)
% 					rmsePlot(:,i)=seqPred(i).Snap(maxSnap).RmseT;		% accessing each of the 50 replication
% 				end
% 				
% 				rmseAvg(:,k)=sum(rmsePlot,2)./(length(seqPred)-(numberCompleteMeasure)); % averaging all replications
% 				
% 			end
% 			
% 			% plotting eig cmparison
% 			lineRmseEig=plot(rmseAvg);
% 			hold on;
% 			
% 		end
% 		
% 		
% 	end
% 	
% 	
% end


