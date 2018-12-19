% Construct Plot Shown in Figure 3
% 12-11-2015, mww


% --- Preliminaries --- 
clear all;
small = 1.0e-10;
big = 1.0e+6;  

% -- File Directories  
outdir = '/Users/mwatson/Dropbox/Hom_Factor/ddisk/Matlab/out/';
figdir = '/Users/mwatson/Dropbox/Hom_Factor/ddisk/Matlab/fig/';
matdir = '/Users/mwatson/Dropbox/Hom_Factor/ddisk/Matlab/mat/';

% --- File Name for output 
outfile_name = [outdir 'Descriptive_statistics_real_variables.out'];
fileID = fopen(outfile_name,'w');

% ----------- Sample Period, Calendars and so forth
[dnobs_m,calvec_m,calds_m] = calendar_make([1959 1],[2014 12],12);  % Monthly Calendar
[dnobs_q,calvec_q,calds_q] = calendar_make([1959 1],[2014 4],4);    % Quarterly Calendar

% -- Load Data
load_data=1;
  % Demeaning Parameters
  i_demean = 1;  % 0 Do Nothing
               % 1 Eliminate low-frequency by local Demeaning
               % 2 Eliminate low-frequency trend by full-sample demeaning;
    
  bw_bw = 100;   % Bi-Weight Parameter for local demeaning
datain_real;        % datain_2 reads in the real-variable dataset

% Calendar, sample, parameters
est_par.smpl_par.calvec = datain.calvec;    % calendar
est_par.smpl_par.nper   = 4;                % number of periods a year

% Factor analysis parameters
% Factor analysis parameters
nfac.unobserved = 1;
nfac.observed = 0;
nfac.total = nfac.unobserved+nfac.observed;
est_par.fac_par.nfac = nfac;
est_par.fac_par.lambda_constraints_est  = 1;      % no constraints on lambda
est_par.fac_par.nt_min                  = 20;     % min number of obs for any series used to est factors
est_par.fac_par.tol                     = 10^-8;  % precision of factor estimation (scaled by by n*t)
est_par.fac_par.nvar_lag                    = 4;      % Number of VAR lags for factors in computign dynamic factor model        

% Factor Parameters
nfac = 1;

% Extract Estimation Sample 
est_data = datain.bpdata(:,datain.bpinclcode==1);
est_namevec = datain.bpnamevec(:,datain.bpinclcode==1);

% Estimate Factor
est_par.smpl_par.nfirst = [1959 3];         % start date
est_par.smpl_par.nlast  = [2014 4];         % end date
lsout_fs= factor_estimation_ls(est_data, est_par);

est_par.smpl_par.nfirst = [1959 3];         % start date
est_par.smpl_par.nlast  = [1983 4];         % end date
lsout_s1= factor_estimation_ls(est_data, est_par);

est_par.smpl_par.nfirst = [1984 1];         % start date
est_par.smpl_par.nlast  = [2014 4];         % end date
lsout_s2= factor_estimation_ls(est_data, est_par);

% Relabel
fac_fs = lsout_fs.fac;
fac_s1 = lsout_s1.fac;
fac_s2 = lsout_s2.fac;

% Normalize split sample values
ii = isnan(fac_s1);
tmp1 = fac_fs(ii==0);
tmp2 = fac_s1(ii==0);
s1 = std(tmp1);
m1 = mean(tmp1);
m2 = mean(tmp2);
s2 = std(tmp2);
tmp1 = tmp1-m1;
tmp2 = ((tmp2-m2)*(s1/s2));
if sign(tmp1(1)) ~= sign(tmp2(1));
    tmp2 = -tmp2;
end;
tmp2 = tmp2+m1;
fac_s1(ii==0)=tmp2;

ii = isnan(fac_s2);
tmp1 = fac_fs(ii==0);
tmp2 = fac_s2(ii==0);
s1 = std(tmp1);
m1 = mean(tmp1);
m2 = mean(tmp2);
s2 = std(tmp2);
tmp1 = tmp1-m1;
tmp2 = ((tmp2-m2)*(s1/s2));
if sign(tmp1(1)) ~= sign(tmp2(1));
    tmp2 = -tmp2;
end;
tmp2 = tmp2+m1;
fac_s2(ii==0)=tmp2;

fac_fs = sum(lagmatrix(fac_fs,0:3),2);
fac_s1 = sum(lagmatrix(fac_s1,0:3),2);
fac_s2 = sum(lagmatrix(fac_s2,0:3),2);

fac_fs = -fac_fs;
fac_s1 = -fac_s1;
fac_s2 = -fac_s2;
    
 plot(datain.calvec,fac_fs,'- k','LineWidth',2);
 hold on;
   plot(datain.calvec,fac_s1,': b','LineWidth',2);
   plot(datain.calvec,fac_s2,'-- r','LineWidth',2);
  hold off;
  %titstr = char(splot.title(i));;
  %title(titstr,'FontSize',18);
  xlim([1955 2020]);
  legend('Full-Sample Estimate','1959-83 Estimate','1984-2014 Estimate');
  legend('Location','NorthWest');
  ax = gca;
  ax.FontSize = 20;
  figname = [figdir 'figure_3'];
  savefig(figname);
  
  e = packr([fac_fs fac_s1]);
  nt = size(e,1);
  emean = mean(e)';   % mean 
  estd = std(e)';     % std 
  ebal = (e - repmat(emean',nt,1))./repmat(estd',nt,1);   % standardized data
  cor_1 = ebal'*ebal/(nt-1); 
  
  e = packr([fac_fs fac_s2]);
  nt = size(e,1);
  emean = mean(e)';   % mean 
  estd = std(e)';     % std 
  ebal = (e - repmat(emean',nt,1))./repmat(estd',nt,1);   % standardized data
  cor_2 = ebal'*ebal/(nt-1); 
  



