% Tabulate Variance Decomps (Dynamic Factors) and R-squared (Static Factors) 
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

str_out = 'tabulate_variance_decomps_';

% -- Load Data
load_data=1;
  % Demeaning Parameters
  i_demean = 1;  % 0 Do Nothing
               % 1 Eliminate low-frequency by local Demeaning
               % 2 Eliminate low-frequency trend by full-sample demeaning;
    
  bw_bw = 100;   % Bi-Weight Parameter for local demeaning
datain_all;      % datain_all reads in the full dataset .. all variables, etc. saved in datain.xx
n_series = size(datain.bpdata,2);

% Factors
% Number of unobserved factors
nfac.unobserved = 8;
% Number of bbserved Factors ..
nfac.observed = 0;
nfac.total = nfac.unobserved+nfac.observed;
est_par.fac_par.nfac = nfac;
est_par.fac_par.w = 1;          % This contains data on the observed factors (scalar if none)


% Sampling parameters
est_par.smpl_par.nfirst = [1959 3];       % start date
est_par.smpl_par.nlast  = [2014 4];       % end date
est_par.smpl_par.calvec = datain.calvec;  % calendar
est_par.smpl_par.nper   = 4;              % number of periods a year

% Factor analysis parameters
est_par.fac_par.nt_min                  = 20;     % min number of obs for any series used to est factors
est_par.lambda.nt_min                   = 40;     % min number of obs for any series used to estimate lamba, irfs, etc.
est_par.fac_par.tol                     = 10^-8;  % precision of factor estimation (scaled by by n*t)

%   -- Get List of variable names
est_namevec = datain.bpnamevec(:,datain.bpinclcode==1);

% Restrictions on factor loadings to identify factors
est_par.fac_par.lambda_constraints_est  = 1;  % no constraints on lambda
est_par.fac_par.lambda_constraints_full = 1;  % no constraints on lambda

% VAR parameters
est_par.var_par.nlag   = 4;    % number of lags
est_par.var_par.iconst = 1;    % include constant
est_par.var_par.icomp  = 1;    % compute companion form of model .. excluding constant

% yit equation parameters
est_par.n_uarlag = 4;  % number of arlags for uniqueness

% Decomposition parameters
decomp_par.varcum = 1;   % use cumulative variance (that is, fraction of variance from factors 1-j, for j = 1, ... 
decomp_par.cancor = 1;   % use cannonical correlations for determining dynamic factors  ... reduced form analysis with dyn. factors ordered by can canonical correlations with data 
decomp_par.hor    = 13;  % horizon for impulse responses and variance decomposition

% Estimate dynamic factor model
fac_est_out = factor_estimation_ls_full(datain.bpdata,datain.bpinclcode,est_par);                          % estimation
irf_vdecomp_out = dynamic_factor_irf_vdecomp(fac_est_out, est_par, decomp_par, datain.bptcodevec);         % irf and variance decomposition


% ------------------------------------------------------------------------

  % Save VDs in CSV File
  % Write Results
  nimp=decomp_par.hor;
  str_out= [outdir str_out num2str(nfac.total) '.csv'];
  fid_out =  fopen(str_out,'w');
  fprintf(fid_out,'Variance decomps for dynamic factors \n');
  fprintf(fid_out,'  Number of static factors: %2i \n',nfac.total);
  fprintf(fid_out,'    Number of VAR lags: %2i \n',est_par.var_par.nlag);
  fprintf(fid_out,'    Number of AR lags for uniqueness: %2i \n',est_par.n_uarlag);
  fprintf(fid_out,'    Each row shows cummulative variance decomps for n_dyn factors\n');
  fprintf(fid_out,'    Columns denote horizons (as in IRFs) \n');
  fprintf(fid_out,'\n\n');
  fprintf(fid_out,'Name,Desc,tcode,n_dyn,,');
  tmp = (0:1:nimp-1);
  prtmat_comma(tmp,fid_out,'%2i','\n');
  for is = 1:n_series;
   tmp = char(datain.bpnamevec(is));
   fprintf(fid_out,[tmp ',']);
   tmp = char(datain.bplabvec_short(is));
   fprintf(fid_out,[tmp ',']);
   tc = datain.bptcodevec(is);
   fprintf(fid_out,'%2i ,',tc);
   fprintf(fid_out,'1,,');
   tmp = reshape(irf_vdecomp_out.vfrac_y_fac_mat(is,1,:),1,nimp); 
   prtmat_comma(tmp,fid_out,'%6.3f','\n');
   for ii = 2:nfac.total;
     fprintf(fid_out,',,,%2i ,,',ii);  
     tmp = reshape(irf_vdecomp_out.vfrac_y_fac_mat(is,ii,:),1,nimp); 
     prtmat_comma(tmp,fid_out,'%6.3f','\n');
   end;
  end;
