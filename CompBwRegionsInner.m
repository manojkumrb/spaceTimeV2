%% Kalman with krigging the error term, including functions for interpolating eigen values for inner
% with aim to plot performance chart according to regions

close all;
clear ;

run 'C:\Users\babu_m\Documents\GitHub\source\initFiles.m';

%% load pre-defined part regions
partRegions     =load('inner_regions.mat');% contains variable nodeIDoutRect
nodeIDoutRectDense   =partRegions.nodeIDoutRect;

devPatterns     =load('devAutoCorrInner_dense.mat');%simAutoCorDevInnerBatchesCombined.mat');%  %hingeDevArSimLoc
devPatterns     =devPatterns.femDevDomain;

% devPatterns(l).FemDevDomain
%% loading from predefined selection
selNodes=load('doorInnerSelNodes.mat');
selNodes=selNodes.selNodes;
Mnp=selNodes.Coord;
iMnp=selNodes.ID;

%% finding key points in each region
keyPointIndex		  = zeros(size(nodeIDoutRectDense,1),1);
keyPointIndex(iMnp)	  = 1;
nodeIDoutRectCoarse= bsxfun(@times,keyPointIndex,nodeIDoutRectDense);


%%  calculating region wise rmse

% reading various predicted deviations
nBasis		=30;%[2,7,30,60];			% set containing number of basis vectors
nBatches	=[1,2];					% set containing number of batches trained on
numberCompleteMeasure=10;
maxSnap		=20;


for k=1:length(nBasis)
	
	fileString=sprintf('allDataInner_eig%i.mat',nBasis(k));
	seqPred=load (fileString);
	seqPred=seqPred.seqPred;
	
	for i=numberCompleteMeasure+1:numberCompleteMeasure+2 %length(seqPred)
		
		devSim=devPatterns(i,:);
		
		bar1=figure;
		bar2=figure;
		
		for ii=1:maxSnap
			
			devPredict=seqPred(i).Snap(ii).Dev;
			
			for iii=1:size(nodeIDoutRectDense,2)
				regionNodeIDs=find(nodeIDoutRectDense(:,iii));
				devPredRegion{iii} = devPredict(regionNodeIDs);
				devSimRegion{iii}  = devSim(regionNodeIDs);
				diff=devPredRegion{iii}-devSimRegion{iii};
				rmsePlot(ii,iii)=sqrt(sum((diff).^2)/size(devPredRegion{iii},2));  % rmse without noisy input
				meanDevRegion(ii,iii)=mean(devSimRegion{iii});

			end
			
			if ii<=size(nodeIDoutRectDense,2)/2
				set(0,'currentFigure',bar1)
				subplot(size(nodeIDoutRectDense,2)/4,2,ii)
				barAxis=bar(rmsePlot(ii,:),0.4);
				% 			barAxis.FaceColor='r';
				hold all
				barAxis2=bar(meanDevRegion(ii,:),0.6);
				barAxis2.FaceAlpha=0.1;
				titlestring= sprintf('Measured region %i in measurement number %i',seqPred(i).Snap(maxSnap).Region(ii),ii);
				ntitle(titlestring,'fontsize',12,'fontname','Times New Roman');
				% 			barAxis2.FaceColor='r';
% 				legend('RMSE','Mean Deviation');
				hold off
			else
				set(0,'currentFigure',bar2)
				indexSubAxis=mod(ii,size(nodeIDoutRectDense,2)/2);
				if indexSubAxis==0
					indexSubAxis=size(nodeIDoutRectDense,2)/2;
				end
				subplot(size(nodeIDoutRectDense,2)/4,2,indexSubAxis)
				barAxis=bar(rmsePlot(ii,:),0.4);
				% 			barAxis.FaceColor='r';
				hold all
				barAxis2=bar(meanDevRegion(ii,:),0.6);
				barAxis2.FaceAlpha=0.1;
				titlestring= sprintf('Measured region %i in measurement number %i',seqPred(i).Snap(maxSnap).Region(ii),ii);
				ntitle(titlestring,'fontsize',12,'fontname','Times New Roman');
				% 			barAxis2.FaceColor='r';
% 				legend('RMSE','Mean Deviation');
				hold off	
				
				
			end
			

		end
		set(0,'currentFigure',bar1)
		stringTitle='Region wise RMSE comparison after each measurement';
		stringXlabel='Region Number';
		stringYlabel='mm';
		figtitle(stringTitle,'fontsize',14,'fontname','Times New Roman')
		bar1=changeAxesLooksBarChart(bar1,stringTitle,stringXlabel,stringYlabel);
		
		set(0,'currentFigure',bar2)
		stringTitle='Region wise RMSE comparison after each measurement';
		stringXlabel='Region Number';
		stringYlabel='mm';
		figtitle(stringTitle,'fontsize',14,'fontname','Times New Roman')
		bar2=changeAxesLooksBarChart(bar2,stringTitle,stringXlabel,stringYlabel);
	end
		
end





