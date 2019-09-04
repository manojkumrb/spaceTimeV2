function fig=plotRmseAndWrong(rmseT,wID,fig)
% plots the time series of state variable and its variance 
% in a uiPanel, grabs data from the figure for the previous timestep
% att	= column vector of state variables
% ptt	= column vector of state variables' variances
% fig	= handle to the figure containing the previous plot, if empty plots
% a new figure.
% 28/05/2018

if nargin <3
	scrsz		= get(groot,'ScreenSize');
	fig			= figure('OuterPosition',[scrsz(3)/6 scrsz(4)/6 scrsz(3)*4/6 scrsz(4)*4/6]);
	
	panhandle	= uipanel('parent',fig,'Position', [0 0 1 .5]);
	
	axWrong		= subplot(2,1,1,'Parent', panhandle,'tag','axWrong', 'NextPlot', 'add');
	plot(axWrong,wID,'*','Tag','wrongID') ;
	axWrong.XTick =1:9;
	axWrong.XLabel.String	='Measurement Number';
	axWrong.YLabel.String	={'Percentage of','correct predictions'};

	axRmse					= subplot(2,1,2,'Parent', panhandle,'Tag','axRmse', 'NextPlot', 'add');
	plot(axRmse,rmseT,'*','tag','rmse') ;
	axRmse.XTick =1:9;
	axRmse.XLabel.String	='Measurement Number';
	axRmse.YLabel.String	='RMSE in mm';
	return
else
	gcf			= fig;
	wrongID		= findobj(gcf,'tag','wrongID');
	rmse		= findobj(gcf,'tag','rmse');
end



oldID		= get(wrongID,'YData'); % gives an cell array (lines X 1)
oldrmse		= get(rmse,'YData');

newID		=[oldID, wID];
newrmse		=[oldrmse, rmseT];



axWrong		=findobj(fig,'Tag','axWrong');
cla(axWrong)

axRmse		=findobj(fig,'Tag','axRmse');
cla(axRmse)
plot(axWrong,newID','tag','wrongID') ;
plot(axRmse,newrmse','tag','rmse') ;
drawnow;


end