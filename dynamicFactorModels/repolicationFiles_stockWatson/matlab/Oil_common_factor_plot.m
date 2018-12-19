% Estimate SDFM.
% IRFS and VDs are for Cholesky Ordered shocks in Factor Model
% Factors can be identified using restrictions on factor loadings
% November 23, mww and po
% 

clear all;
small = 1.0e-10;
big = 1.0e+6;
rng(63287545);

% -- File Directories  
outdir = '/Users/mwatson/Dropbox/Hom_Factor/ddisk/Matlab/out/';
figdir = '/Users/mwatson/Dropbox/Hom_Factor/ddisk/Matlab/fig/';
matdir = '/Users/mwatson/Dropbox/Hom_Factor/ddisk/Matlab/mat/'; 

% ---- Identifier for output files and figures ----
% str_ouput = 'timing_restrictions_oil_price_2_fs';
% str_ouput = 'timing_restrictions_oil_price_2_ps1';
str_ouput = 'temp';

% ----------- Sample Period, Calendars and so forth
[dnobs_m,calvec_m,calds_m] = calendar_make([1959 1],[2014 12],12);
[dnobs_q,calvec_q,calds_q] = calendar_make([1959 1],[2014 4],4);

% -- Load Data
load_data=1;
  % Demeaning Parameters
  i_demean = 1;  % 0 Do Nothing
               % 1 Eliminate low-frequency by local Demeaning
               % 2 Eliminate low-frequency trend by full-sample demeaning;
    
  bw_bw = 100;   % Bi-Weight Parameter for local demeaning
datain_all;
n_series = size(datain.bpdata,2);

% Factors
% Number of unobserved factors
nfac.unobserved = 8;
% .. Observed Factors ..
nfac.observed = 0;
nfac.total = nfac.unobserved+nfac.observed;
est_par.fac_par.nfac = nfac;
est_par.fac_par.w = 1;         % This contains data on the observed factors (scalar if none)


% Sampling parameters
est_par.smpl_par.nfirst = [1985 1];       % start date
est_par.smpl_par.nlast  = [2014 4];       % end date
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
str_var = {'WPU0561', 'MCOILWTICO', 'MCOILBRENTEU', 'RAC_IMP'};
R = eye(nfac.total);
r = [1; zeros(nfac.total-1,1)];
lambda_constraints_estdata = lambda_construct(str_var, est_namevec, R,r);
lambda_constraints_bpdata = lambda_construct(str_var, datain.bpnamevec, R,r);% other oil variables that load on only first factor, but without unit-scale restriction
est_par.fac_par.lambda_constraints_est = lambda_constraints_estdata;
est_par.fac_par.lambda_constraints_full = lambda_constraints_bpdata;

% VAR parameters
est_par.var_par.nlag   = 4;  % number of lags
est_par.var_par.iconst = 1;  % include constant
est_par.var_par.icomp  = 1;  % compute Companion form of model .. excluding constant

% yit equation parameters
est_par.n_uarlag = 4;  % Number of arlags for uniqueness

% Decomposition parameters
decomp_par.varcum = 0;   % don't use cumulative variance
decomp_par.cancor = 0;   % don't use cannonical correlations
decomp_par.hor    = 17;  % horizon for impulse responses and variance decomposition


% -- Construct Estimates of Factors --- 
fac_est_out = factor_estimation_ls_full(datain.bpdata,datain.bpinclcode,est_par);                  % estimation

str_var = {'WPU0561', 'MCOILWTICO', 'MCOILBRENTEU', 'RAC_IMP'};
ustr = char(str_var);
j_oilprice = colnumber(ustr,datain.bpnamevec);
% -- Get Series and Fitted Values
oilprice_data = datain.bpdata(:,j_oilprice);
lam_oil = fac_est_out.lam_mat(j_oilprice,:);
oilp_fv = fac_est_out.fac*lam_oil';
oilp_fv = oilp_fv(:,1);
op_4q = NaN(size(oilprice_data));
opf_4q = NaN(size(oilp_fv));
tmp = [oilprice_data oilp_fv];
op_4q(4:end,:) = 100*(oilprice_data(1:end-3,:)+oilprice_data(2:end-2,:)+oilprice_data(3:end-1,:)+oilprice_data(4:end,:));
opf_4q(4:end,:) = 100*(oilp_fv(1:end-3,:)+oilp_fv(2:end-2,:)+oilp_fv(3:end-1,:)+oilp_fv(4:end,:));


% ---- Figure 7b ---
plot(datain.calvec,400*oilprice_data,'LineWidth',2);
hold on;
  plot(datain.calvec,400*oilp_fv,': k','LineWidth',4);
hold off;
legend('Oil Price: PPI','Oil Price: WTI','Oil Price: Brent','Oil Price: RAC','Common Component');
legend('Location','NorthEast');
xlabel('Date');
ylabel('PAAR');
xlim([1985 2015]);
ax = gca;
ax.FontSize = 20;
figname = [figdir 'figure_7b'];
savefig(figname); 

% ---- Figure 7a -----
str='MCOILBRENTEU';    % Brent Oil Price 
str=upper(str);
j = colnumber(str,datain.bpnamevec);
poil_brent = datain.bpdata_raw(:,j);
poil_brent = 100*poil_brent;
plot(datain.calvec,poil_brent,'LineWidth',2);
xlabel('Date');
ylabel('Dollars ($2009) per Barel');
xlim([1985 2015]);
ax = gca;
ax.FontSize = 20;
figname = [figdir 'figure_7a'];
savefig(figname); 

