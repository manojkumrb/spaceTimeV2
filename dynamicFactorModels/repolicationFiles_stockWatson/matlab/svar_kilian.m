% Estimate SDFM.
% IRFS and VDs are for Cholesky Ordered shocks in Factor Model
% Factors can be identified using restrictions on factor loadings
% February 28, 2016, mww and po
% 
% -- THIS VERSION:  Killian identification.  Uses all observed factors (SVAR).

clear all;
small = 1.0e-10;
big = 1.0e+6;
rng(63287545);

% -- File Directories  
outdir = '/Users/mwatson/Dropbox/Hom_Factor/ddisk/Matlab/out/';
figdir = '/Users/mwatson/Dropbox/Hom_Factor/ddisk/Matlab/fig/';
matdir = '/Users/mwatson/Dropbox/Hom_Factor/ddisk/Matlab/mat/'; 

% ---- Identifier for output files and figures ----
str_ouput = 'SVAR_kilian';
% -- Number of Identified Shocks 
n_ident_shocks = 3;

% ----------- Sample Period, Calendars and so forth
[dnobs_m,calvec_m,calds_m] = calendar_make([1959 1],[2014 12],12);
[dnobs_q,calvec_q,calds_q] = calendar_make([1959 1],[2014 4],4);

% -- Load Data
load_data=0;
  % Demeaning Parameters
  i_demean = 1;  % 0 Do Nothing
               % 1 Eliminate low-frequency by local Demeaning
               % 2 Eliminate low-frequency trend by full-sample demeaning;
    
  bw_bw = 100;   % Bi-Weight Parameter for local demeaning
datain_all;
n_series = size(datain.bpdata,2);

% Factors
% Number of unobserved factors
nfac.unobserved = 0;
% .. Observed Factors ..
str_var = {'OILPROD_SA';'GLOBAL_ACT';'WPU0561';'GDPC96';'PAYEMS';'PCECTPI';'FEDFUNDS';'TWEXMMTH'};
observed_factor = observed_factor_setup(str_var,datain);
nfac.observed = size(observed_factor.fac,2);
est_par.fac_par.w = observed_factor.fac;         % This contains data on the observed factors
nfac.total = nfac.unobserved+nfac.observed;
est_par.fac_par.nfac = nfac;

% Sampling parameters
est_par.smpl_par.nfirst = [1985 1];       % start date
est_par.smpl_par.nlast  = [2014 3];       % end date
est_par.smpl_par.calvec = datain.calvec;  % calendar
est_par.smpl_par.nper   = 4;              % number of periods a year

% Factor analysis parameters
est_par.fac_par.nt_min                  = 20;     % min number of obs for any series used to est factors
est_par.lambda.nt_min                   = 40;     % min number of obs for any series used to estimate lamba, irfs, etc.
est_par.fac_par.tol                     = 10^-8;  % precision of factor estimation (scaled by by n*t)

% Restrictions on Factor Loadings used to identify factors
%   * note order of factors: observed factors and then unobserved factors
%   -- Get List of variable names
est_namevec = datain.bpnamevec(:,datain.bpinclcode==1);

% Identification here is for oil price variables
% First set of variables -- unit coefficients for first factor loading .. others are zeros
%  Note, this identifies first factor as oil price factor with scale determined by these observatbles
%
%
% -- Oil Production is first shock
str_var = {'OILPROD_SA'};
R = eye(nfac.total);
r = [1; zeros(nfac.total-1,1)];
lambda_constraints_estdata_1 = lambda_construct(str_var, est_namevec, R,r);
lambda_constraints_bpdata_1 = lambda_construct(str_var, datain.bpnamevec, R,r);

% -- Global Activity is second shock
str_var = {'GLOBAL_ACT'};
R = eye(nfac.total);
r = [0; 1; zeros(nfac.total-2,1)];
lambda_constraints_estdata_2 = lambda_construct(str_var, est_namevec, R,r);
lambda_constraints_bpdata_2 = lambda_construct(str_var, datain.bpnamevec, R,r);

% -- Oil Price is third shock
str_var = {'WPU0561'};
R = eye(nfac.total);
r = [0; 0; 1; zeros(nfac.total-3,1)];
lambda_constraints_estdata_3 = lambda_construct(str_var, est_namevec, R,r);
lambda_constraints_bpdata_3 = lambda_construct(str_var, datain.bpnamevec, R,r);

% combine constraints
est_par.fac_par.lambda_constraints_est = [lambda_constraints_estdata_1; lambda_constraints_estdata_2; lambda_constraints_estdata_3];
est_par.fac_par.lambda_constraints_full = [lambda_constraints_bpdata_1; lambda_constraints_bpdata_2; lambda_constraints_bpdata_3];


% VAR parameters
est_par.var_par.nlag   = 4;  % number of lags
est_par.var_par.iconst = 1;  % include constant
est_par.var_par.icomp  = 1;  % compute Companion form of model .. excluding constant

% yit equation parameters
est_par.n_uarlag = 4;  % Number of arlags for uniqueness

% Decomposition parameters
decomp_par.varcum = 0;   % don't use cumulative variance -- here show VDs, shock-by-shock (here one shock)
decomp_par.cancor = 0;   % don't use cannonical correlations -- identification is via "named factor" etc.
decomp_par.hor    = 17;  % horizon for impulse responses and variance decomposition


% -- Construct Estimates of Factors --- 
fac_est_out = factor_estimation_ls_full(datain.bpdata,datain.bpinclcode,est_par);                  % estimation
irf_vdecomp_out = dynamic_factor_irf_vdecomp(fac_est_out,est_par,decomp_par,datain.bptcodevec);      % variance decomposition

% -- Compute Standard Errors for IRFs and VDs using parametric bootstrap simulations
n_rep = 500;   % Number of bootstap simulations for computing SEs
se_irf_vdecomp_out = se_dynamic_factor_irf_vdecomp(datain,fac_est_out,est_par,decomp_par,n_rep); 

for ishock = 1:n_ident_shocks;

  % Save IRFs and VDs to CSV File
  % Write Results
  nimp=decomp_par.hor;
  str_out= [outdir str_ouput '_irf_vd_shock' num2str(ishock) '.csv'];
  fid_out =  fopen(str_out,'w');
  fprintf(fid_out,'Variance decomps for dynamic factors \n');
  fprintf(fid_out,'  Number of unobserved static = number of dynamic factors: %2i \n',est_par.fac_par.nfac.unobserved);
  fprintf(fid_out,'  Number of observed static = number of dynamic factors: %2i \n',est_par.fac_par.nfac.observed);
  fprintf(fid_out,'  Identified Shock Number: %2i \n',ishock);
  if est_par.fac_par.nfac.observed > 0;
    fprintf(fid_out,'     Observed static factors:\n');
    for i = 1:est_par.fac_par.nfac.observed;
      tmp = char(observed_factor.names(i));
      fprintf(fid_out,['          ' tmp '\n']);
    end;
  end;
  fprintf(fid_out,'     Sample Period: %4i:Q%1i',est_par.smpl_par.nfirst);
  fprintf(fid_out,'-%4i:Q%1i\n',est_par.smpl_par.nlast);
  fprintf(fid_out,'    Number of VAR lags: %2i \n',est_par.var_par.nlag);
  fprintf(fid_out,'    Number of AR lags for uniqueness: %2i \n',est_par.n_uarlag);
  fprintf(fid_out,'    Number of bootstrap simulations for SE calculation: %4i \n',n_rep);
  fprintf(fid_out,'    Columns denote horizons (as in IRFs) \n');
  fprintf(fid_out,'\n\n');
  fprintf(fid_out,'Name,Desc,,');
  tmp = (0:1:nimp-1);
  prtmat_comma(tmp,fid_out,'%2i','\n');
  for is = 1:n_series;
   tmp = char(datain.bpnamevec(is));
   fprintf(fid_out,[tmp ',']);
   tmp = char(datain.bplabvec_short(is));
   fprintf(fid_out,[tmp ',IRF Estimate,']);
   tmp = reshape(irf_vdecomp_out.imp_y_fac_mat_scl(is,ishock,:),1,nimp); 
   prtmat_comma(tmp,fid_out,'%0.4e','\n');
   tmp = reshape(se_irf_vdecomp_out.mean_imp_y_fac_mat_scl(is,ishock,:),1,nimp); 
   fprintf(fid_out,',,IRF bootstrap mean,');
   prtmat_comma(tmp,fid_out,'%0.4e','\n');
   tmp = reshape(se_irf_vdecomp_out.se_imp_y_fac_mat_scl(is,ishock,:),1,nimp); 
   fprintf(fid_out,',,IRF bootstrap std. deviation,');
   prtmat_comma(tmp,fid_out,'%0.4e','\n');
   tmp = reshape(irf_vdecomp_out.vfrac_y_comp_mat(is,ishock,:),1,nimp); 
   fprintf(fid_out,',,VD Estimate,');
   prtmat_comma(tmp,fid_out,'%6.3f','\n');
   tmp = reshape(se_irf_vdecomp_out.mean_vfrac_y_comp_mat(is,ishock,:),1,nimp); 
   fprintf(fid_out,',,VD bootstrap mean,');
   prtmat_comma(tmp,fid_out,'%6.3f','\n');
   tmp = reshape(se_irf_vdecomp_out.se_vfrac_y_comp_mat(is,ishock,:),1,nimp); 
   fprintf(fid_out,',,VD bootstrap std. deviation,');
   prtmat_comma(tmp,fid_out,'%6.3f','\n');  
  end;

end;

 % Variables to Save for future figures and tables
  rslt.datain = datain;
  rslt.irf_vdecomp_out = irf_vdecomp_out;
  rslt.se_irf_vdecomp_out = se_irf_vdecomp_out;
  str_tmp = [matdir 'rslt_' str_ouput]; save(str_tmp,'rslt');

  