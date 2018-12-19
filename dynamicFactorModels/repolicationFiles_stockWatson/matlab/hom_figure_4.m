% Construct Plot Shown in Figure 4
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

% Sampling parameters
est_par.smpl_par.nfirst = [1959 3];         % start date
est_par.smpl_par.nlast  = [2014 4];         % end date
est_par.smpl_par.calvec = datain.calvec;    % calendar
est_par.smpl_par.nper   = 4;                % number of periods a year

% Factor analysis parameters
est_par.fac_par.lambda_constraints_est  = 1;      % no constraints on lambda
est_par.fac_par.nt_min                  = 20;     % min number of obs for any series used to est factors
est_par.lambda.nt_min                   = 40;     % min number of obs for any series used to estimate lamba, irfs, etc.
est_par.fac_par.tol                     = 10^-8;  % precision of factor estimation (scaled by by n*t)
est_par.fac_par.nvar_lag                = 4;      % Number of VAR lags for factors in computign dynamic factor model        

est_par.fac_par.nt_min                  = 20;     % min number of obs for any series used to est factors
est_par.fac_par.tol                     = 10^-8;  % precision of factor estimation (scaled by by n*t)

% Restrictions on factor loadings to identify factors
est_par.fac_par.lambda_constraints_est  = 1;  % no constraints on lambda
est_par.fac_par.lambda_constraints_full = 1;  % no constraints on lambda

% VAR parameters
est_par.var_par.nlag   = 4;    % number of lags
est_par.var_par.iconst = 1;    % include constant
est_par.var_par.icomp  = 1;    % compute companion form of model .. excluding constant

% yit equation parameters
est_par.n_uarlag = 4;  % number of arlags for uniqueness

% Estimate dynamic factor model
nfac.unobserved = 1;
nfac.observed = 0;
nfac.total = nfac.unobserved+nfac.observed;
est_par.fac_par.nfac = nfac;
fac_est_out_1 = factor_estimation_ls_full(datain.bpdata,datain.bpinclcode,est_par);           % estimation
nfac.unobserved = 3;
nfac.observed = 0;
nfac.total = nfac.unobserved+nfac.observed;
est_par.fac_par.nfac = nfac;
fac_est_out_3 = factor_estimation_ls_full(datain.bpdata,datain.bpinclcode,est_par);           % estimation
nfac.unobserved = 5;
nfac.observed = 0;
nfac.total = nfac.unobserved+nfac.observed;
est_par.fac_par.nfac = nfac;
fac_est_out_5 = factor_estimation_ls_full(datain.bpdata,datain.bpinclcode,est_par);           % estimation

j = colnumber('GDPC96',datain.bpnamevec);
y = datain.bpdata(:,j);
y_ma4 = sum(lagmatrix(y,[0 1 2 3]),2);

yf_1 = fac_est_out_1.fac*fac_est_out_1.lam_mat(j,:)';
yf_1_ma4 = sum(lagmatrix(yf_1,[0 1 2 3]),2);   

yf_3 = fac_est_out_3.fac*fac_est_out_3.lam_mat(j,:)';
yf_3_ma4 = sum(lagmatrix(yf_3,[0 1 2 3]),2);

yf_5 = fac_est_out_5.fac*fac_est_out_5.lam_mat(j,:)';
yf_5_ma4 = sum(lagmatrix(yf_5,[0 1 2 3]),2);   
    
plot(datain.calvec,y_ma4,'- k','LineWidth',2);
hold on;
  plot(datain.calvec,yf_1_ma4,'-- b','LineWidth',2);
  plot(datain.calvec,yf_3_ma4,': r','LineWidth',2);
  plot(datain.calvec,yf_5_ma4,'-. c','LineWidth',2);
hold off;
xlim([1955 2020]);
legend('GDP Growth','1-Factor Common Component','3-Factor Common Component','5-Factor Common Component');
legend('Location','NorthWest');
ax = gca;
ax.FontSize = 20;
figname = [figdir 'figure_4'];
savefig(figname);

% Compute R2 Values
tmp = packr([y_ma4 yf_1_ma4 yf_3_ma4 yf_5_ma4]);
y = tmp(:,1);
x = tmp(:,2);
y = y - mean(y);
x = x - mean(x);
b = (y'*x)/(x'*x);
e = y - x*b;
sse = e'*e;
ssy = y'*y;
r2_1 = 1-(sse/ssy);

x = tmp(:,3);
y = y - mean(y);
x = x - mean(x);
b = (y'*x)/(x'*x);
e = y - x*b;
sse = e'*e;
ssy = y'*y;
r2_3 = 1-(sse/ssy);

x = tmp(:,4);
y = y - mean(y);
x = x - mean(x);
b = (y'*x)/(x'*x);
e = y - x*b;
sse = e'*e;
ssy = y'*y;
r2_5 = 1-(sse/ssy);
