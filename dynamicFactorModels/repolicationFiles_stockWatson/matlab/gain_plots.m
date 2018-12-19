% Gain plots
%

% --- Preliminaries --- 
clear all;
small = 1.0e-10;
big = 1.0e+6;  

% -- File Directories  
outdir = '/Users/mwatson/Dropbox/Hom_Factor/ddisk/Matlab/out/';
figdir = '/Users/mwatson/Dropbox/Hom_Factor/ddisk/Matlab/fig/';
matdir = '/Users/mwatson/Dropbox/Hom_Factor/ddisk/Matlab/mat/';

% --- File Name for output 
outfile_name = [outdir 'Descriptive_statistics_all_variables.out'];
fileID = fopen(outfile_name,'w');

% Frequencies;
wvec = linspace(0,pi,500)';

% Bi-weight paramter   
bw_bw = 100;   % Bi-Weight Parameter for local demeaning
trend = (0:1:100)';
dt = trend/bw_bw;
bw_weight = (15/16)*((1-dt.^2).^2);   % Bi-Weight 
bw_weight = bw_weight.*(abs(dt) < 1);
tmp = flipud(bw_weight(2:end));
bw_weight = [tmp ; bw_weight];
bw_weight = bw_weight/sum(bw_weight);
tmp = flipud(trend(2:end));
trend = [-tmp ; trend];

% Compute Gain
bw_gain = NaN(size(wvec,1),1);
for i = 1:size(wvec,1);
  bw_gain(i) = gain(bw_weight,wvec(i));
end;

% Flat over + and - 40 quarters
ma40_weight = abs(trend) <= 40;
ma40_weight = ma40_weight/sum(ma40_weight);
% Compute Gain
ma40_gain = NaN(size(wvec,1),1);
for i = 1:size(wvec,1);
  ma40_gain(i) = gain(ma40_weight,wvec(i));
end;

% Band pass, 100 terms on each side, 200 quarter cutoff;
nper = 200;
ombar = 2*pi/nper;
t1 = (1:1:100)';
tmp0 = ombar/pi;
tmp1 = (1./(pi*t1)).*sin(t1*ombar);
bp_weight = [flipud(tmp1);tmp0;tmp1];
bp_weight = bp_weight/sum(bp_weight);
% Compute Gain
bp_gain = NaN(size(wvec,1),1);
for i = 1:size(wvec,1);
  bp_gain(i) = gain(bp_weight,wvec(i));
end;

% HP
hp_weight = dlmread('hpfilter_trend.asc');
% Compute Gain
hp_gain = NaN(size(wvec,1),1);
for i = 1:size(wvec,1);
  hp_gain(i) = gain(hp_weight,wvec(i));
end;

% Filtered Weights and Series
bwc_weight = -bw_weight; bwc_weight(101) = bwc_weight(101)+1;
hpc_weight = -hp_weight; hpc_weight(101) = hpc_weight(101)+1;
bpc_weight = -bp_weight; bpc_weight(101) = bpc_weight(101)+1;
ma40c_weight = -ma40_weight; ma40c_weight(101) = ma40c_weight(101)+1;

bwc_gain = ones(500,1)-bw_gain;
bpc_gain = ones(500,1)-bp_gain;
hpc_gain = ones(500,1)-hp_gain;
ma40c_gain = ones(500,1)-ma40_gain;

% Plot Trend Filter Weights
plot(trend,bw_weight,'- k','LineWidth',2);
hold on;
 plot(trend,ma40_weight,'-- b','LineWidth',2);
 plot(trend,bp_weight,': r','LineWidth',2);
 plot(trend,hp_weight,'-. g','LineWidth',2);
hold off;
legend('Biweight','MA(40)','Bandpass','Hodrick-Prescott');
legend('Location','NorthWest');
xlabel('lag');
ax = gca;
ax.FontSize = 20;
figname = [figdir 'figure_2a'];
savefig(figname);

plot(wvec,bw_gain,'- k','LineWidth',2);
hold on;
 plot(wvec,ma40_gain,'-- b','LineWidth',2);
 plot(wvec,bp_gain,': r','LineWidth',2);
 plot(wvec,hp_gain,'-. g','LineWidth',2);
hold off;
legend('Biweight','MA(40)','Bandpass','Hodrick-Prescott');
legend('Location','NorthEast');
xlabel('frequency');
ax = gca;
ax.FontSize = 20;
figname = [figdir 'figure_2b'];
savefig(figname);
