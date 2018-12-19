% Tabulate R-squared for statitic factors for all series.  Save as CSV file 
%

clear all;
small = 1.0e-10;
big = 1.0e+6;

% -- File Directories  
outdir = '/Users/mwatson/Dropbox/Hom_Factor/ddisk/Matlab/out/';
figdir = '/Users/mwatson/Dropbox/Hom_Factor/ddisk/Matlab/fig/';
matdir = '/Users/mwatson/Dropbox/Hom_Factor/ddisk/Matlab/mat/';

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
datain_all;      % datain_all reads in the full dataset .. all variables, etc. saved in datain.xx


% Factor Parameters
nfac_max = 10;

% Sampling parameters
est_par.smpl_par.nfirst = [1959 3];       % start date
est_par.smpl_par.nlast  = [2014 4];       % end date
est_par.smpl_par.calvec = datain.calvec;  % calendar
est_par.smpl_par.nper   = 4;              % number of periods a year

% Factor analysis parameters
est_par.fac_par.nt_min                  = 20;     % min number of obs for any series used to est factors
est_par.lambda.nt_min                   = 40;     % min number of obs for any series used to estimate lamba, irfs, etc.
est_par.fac_par.tol                     = 10^-8;  % precision of factor estimation (scaled by by n*t)

% Restrictions on factor loadings to identify factors
est_par.fac_par.lambda_constraints_est  = 1;  % no constraints on lambda
est_par.fac_par.lambda_constraints_full = 1;  % no constraints on lambda

% VAR parameters for factors
est_par.var_par.nlag   = 4;    % number of lags
est_par.var_par.iconst = 1;    % include constant
est_par.var_par.icomp  = 1;    % compute companion form of model .. excluding constant

% yit equation parameters
est_par.n_uarlag = 4;  % number of arlags for uniqueness

% Matrices for storing results
n_series = size(datain.bpdata,2);
r2_mat = NaN(n_series,nfac_max);

for nfac = 1:nfac_max;
 % Estimate dynamic factor model
 est_par.fac_par.nfac.unobserved = nfac;
 est_par.fac_par.nfac.observed = 0;
 est_par.fac_par.nfac.total = nfac;
 fac_est_out = factor_estimation_ls_full(datain.bpdata,datain.bpinclcode,est_par);                  % estimation
 r2_mat(:,nfac) = fac_est_out.r2;
end;

% Save RSquared Values to CSV File
str_out= [outdir 'tabulate_rsquared_static_factors.csv'];
fid_out =  fopen(str_out,'w');
fprintf(fid_out,'Name,Desc,,');
tmp = (1:1:nfac_max);
prtmat_comma(tmp,fid_out,'%2i','\n');
for is = 1:n_series;
   tmp = char(datain.bpnamevec(is));
   fprintf(fid_out,[tmp ',']);
   tmp = char(datain.bplabvec_short(is));
   fprintf(fid_out,[tmp ',,']); 
   prtmat_comma(r2_mat(is,:),fid_out,'%6.3f','\n');
end;


