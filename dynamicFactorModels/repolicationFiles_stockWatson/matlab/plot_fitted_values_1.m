% Plot selected Fitted Values 
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
datain_real;        % datain_2 reads in the real-variable dataset



% PARAMETERS

% Number of unobserved factors
nfac.unobserved = 1;
% .. Observed Factors ..
nfac.observed = 0;
nfac.total = nfac.unobserved+nfac.observed;
est_par.fac_par.nfac = nfac;

% Sampling parameters
est_par.smpl_par.nfirst = [1959 3];         % start date
est_par.smpl_par.nlast  = [2014 4];         % end date
est_par.smpl_par.calvec = datain.calvec;    % calendar
est_par.smpl_par.nper   = 4;                % number of periods a year


% Restrictions on factor loadings to identify factors
est_par.fac_par.lambda_constraints_est  = 1;  % no constraints on lambda
est_par.fac_par.lambda_constraints_full = 1;  % no constraints on lambda

% VAR parameters
est_par.var_par.nlag   = 4;    % number of lags
est_par.var_par.iconst = 1;    % include constant
est_par.var_par.icomp  = 1;    % compute companion form of model .. excluding constant

est_par.fac_par.lambda_constraints_est  = 1;      % no constraints on lambda
est_par.fac_par.nt_min                  = 20;     % min number of obs for any series used to est factors
est_par.lambda.nt_min                   = 40;     % min number of obs for any series used to estimate lamba, irfs, etc.
est_par.fac_par.tol                     = 10^-8;  % precision of factor estimation (scaled by by n*t)
est_par.fac_par.nvar_lag                = 4;      % Number of VAR lags for factors in computign dynamic factor model        



% yit equation parameters
est_par.n_uarlag = 4;  % number of arlags for uniqueness

% Estimate dynamic factor model
fac_est_out = factor_estimation_ls_full(datain.bpdata,datain.bpinclcode,est_par);           % estimation

% List series to plot
splot.name = {...
    'GDPC96' ...
    'INDPRO' ...
    'PAYEMS' ...
    'A0M057' ...
    };
% Figure Titles
splot.title = {...
    'GDP' ...
    'Industrial Production' ...
    'Nonfarm Employment' ...
    'Manufacturing and Trade Sales' ...
    };
% Construct Plot
for i = 1:size(splot.name,2);
    a1 = char(splot.name(i));
    j = colnumber(a1,datain.bpnamevec);
    y = datain.bpdata(:,j);
    yf = fac_est_out.lam_mat(j,:)*fac_est_out.fac;
    y4 = sum(lagmatrix(y,[0 1 2 3]),2);
    yf4 = sum(lagmatrix(yf,[0 1 2 3]),2);
    y4 = 100*y4;
    yf4 = 100*yf4;
    plot(datain.calvec,y4,'- k','LineWidth',2);
    hold on;
    plot(datain.calvec,yf4,'-- r','LineWidth',2);
    hold off;
    titstr = char(splot.title(i));;
    title(titstr,'FontSize',18);
    xlim([1955 2020]);
    if i == 1;
      legend('Series','Common Component');
      legend('Location','NorthEast');
    end;
    ax = gca;
    ax.FontSize = 20;
    figname = [figdir 'fitted_value_1_' num2str(i)];
    savefig(figname);
end;


