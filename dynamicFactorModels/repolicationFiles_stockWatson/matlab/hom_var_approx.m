% HOM_Table_4 .. compute canonical correlations between factors and a set of variables
% 12-112015, mww


% --- Preliminaries --- 
clear all;
small = 1.0e-10;
big = 1.0e+6;  

% -- File Directories  
outdir = '/Users/mwatson/Dropbox/Hom_Factor/ddisk/Matlab/out/';
figdir = '/Users/mwatson/Dropbox/Hom_Factor/ddisk/Matlab/fig/';
matdir = '/Users/mwatson/Dropbox/Hom_Factor/ddisk/Matlab/mat/';

% --- File Name for output 
outfile_name = [outdir 'hom_var_approx.out'];
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
datain_all;        % datain_1 reads in the full dataset

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

% VAR parameters
est_par.var_par.nlag   = 4;    % number of lags
est_par.var_par.iconst = 1;    % include constant
est_par.var_par.icomp  = 1;    % compute companion form of model .. excluding constant

% yit equation parameters
est_par.n_uarlag = 4;  % number of arlags for uniqueness

% Factor Parameters
% Number of unobserved factors
nfac.unobserved = 8;
% .. Observed Factors ..
nfac.observed = 0;
nfac.total = nfac.unobserved+nfac.observed;
est_par.fac_par.nfac = nfac;

% Extract Estimation Sample 
est_data = datain.bpdata(:,datain.bpinclcode==1);
est_namevec = datain.bpnamevec(:,datain.bpinclcode==1);

% Compute Various Statistics for estimating the number of static and dynamic factors
fac_est_out = factor_estimation_ls_full(datain.bpdata,datain.bpinclcode,est_par);                  % estimation
% Rename
fac = fac_est_out.fac;   % factors
fac_res = fac_est_out.varout.resid;


fprintf(fileID,'Canonical Correlations between VAR variables and Factors \n\n');

% ----------------- VAR A --------------
% Variables to Include in VAR
splot.name = {...
    'GDPC96' ...
    'PAYEMS' ...
    'PCECTPI' ...
    'FEDFUNDS' ...
    };
% Get Series
for i = 1:size(splot.name,2);
  a1 = char(splot.name(i));
  j = colnumber(a1,datain.bpnamevec);
  y = datain.bpdata(:,j);
  if i == 1;
      xvar = y;
  else;
      xvar = [xvar y];
  end;
end;
% Compute VAR
varout_xvar = varest(xvar,est_par.var_par,est_par.smpl_par);
xvar_res = varout_xvar.resid;
% Levels of Facots and Xvariables
x = xvar;
y = fac;
ii = sum(isnan(x),2)+sum(isnan(y),2);
x = x(ii==0,:);
y = y(ii==0,:);
[~,~,r_lev] = canoncorr(x,y);   % NOTE -- Matlab demeans the series 

x = xvar_res;
y = fac_res;
ii = sum(isnan(x),2)+sum(isnan(y),2);
x = x(ii==0,:);
y = y(ii==0,:);
[~,~,r_res] = canoncorr(x,y);   % NOTE -- Matlab demeans the series 

fprintf(fileID,'  Variables included in VAR: \n');
for i = 1:size(splot.name,2);
   a1 = char(splot.name(i));
   fprintf(fileID,['     ' a1 '\n']);  
end;
fprintf(fileID,'  Sample Period: %4i:Q%1i',est_par.smpl_par.nfirst);
fprintf(fileID,'-%4i:Q%1i\n',est_par.smpl_par.nlast);
fprintf(fileID,'  Number of Factors: %3i \n',est_par.fac_par.nfac.unobserved);
fprintf(fileID,'    Canonical Correlations for Levels of X and Factor \n      ');
prtmat_comma(r_lev,fileID,'%5.2f','\n');
fprintf(fileID,'    Canonical Correlations for VAR Residuals of X and Factor \n      ');
prtmat_comma(r_res,fileID,'%5.2f','\n');


% ----------------- VAR B --------------
fprintf(fileID,'\n\n\n\n');
% Variables to Include in VAR
splot.name = {...
    'GDPC96' ...
    'PAYEMS' ...
    'PCECTPI' ...
    'FEDFUNDS' ...
    'NAPMPRI'...
    'WPU0561' ... 
    'CP90_TBILL' ...
    'GS10_TB3M' ...
    };
% Get Series
for i = 1:size(splot.name,2);
  a1 = char(splot.name(i));
  j = colnumber(a1,datain.bpnamevec);
  y = datain.bpdata(:,j);
  if i == 1;
      xvar = y;
  else;
      xvar = [xvar y];
  end;
end;
% Compute VAR
varout_xvar = varest(xvar,est_par.var_par,est_par.smpl_par);
xvar_res = varout_xvar.resid;
% Levels of Facots and Xvariables
x = xvar;
y = fac;
ii = sum(isnan(x),2)+sum(isnan(y),2);
x = x(ii==0,:);
y = y(ii==0,:);
[~,~,r_lev] = canoncorr(x,y);   % NOTE -- Matlab demeans the series 

x = xvar_res;
y = fac_res;
ii = sum(isnan(x),2)+sum(isnan(y),2);
x = x(ii==0,:);
y = y(ii==0,:);
[~,~,r_res] = canoncorr(x,y);   % NOTE -- Matlab demeans the series 

fprintf(fileID,'  Variables included in VAR: \n');
for i = 1:size(splot.name,2);
   a1 = char(splot.name(i));
   fprintf(fileID,['     ' a1 '\n']);  
end;
fprintf(fileID,'  Sample Period: %4i:Q%1i',est_par.smpl_par.nfirst);
fprintf(fileID,'-%4i:Q%1i\n',est_par.smpl_par.nlast);
fprintf(fileID,'  Number of Factors: %3i \n',est_par.fac_par.nfac.unobserved);
fprintf(fileID,'    Canonical Correlations for Levels of X and Factor \n      ');
prtmat_comma(r_lev,fileID,'%5.2f','\n');
fprintf(fileID,'    Canonical Correlations for VAR Residuals of X and Factor \n      ');
prtmat_comma(r_res,fileID,'%5.2f','\n');

% ----------------- VAR C (oil) --------------
fprintf(fileID,'\n\n\n\n');
% Variables to Include in VAR
{'OILPROD_SA';'GLOBAL_ACT';'WPU0561';'GDPC96';'PAYEMS';'PCECTPI';'FEDFUNDS';'TWEXMMTH'}
splot.name = {...
    'OILPROD_SA' ...
    'GLOBAL_ACT' ...
    'WPU0561' ...
    'GDPC96' ...
    'PAYEMS' ...
    'PCECTPI' ...
    'FEDFUNDS' ...
    'TWEXMMTH' ...
    };
% Get Series
for i = 1:size(splot.name,2);
  a1 = char(splot.name(i));
  j = colnumber(a1,datain.bpnamevec);
  y = datain.bpdata(:,j);
  if i == 1;
      xvar = y;
  else;
      xvar = [xvar y];
  end;
end;
% Compute VAR
varout_xvar = varest(xvar,est_par.var_par,est_par.smpl_par);
xvar_res = varout_xvar.resid;
% Levels of Facots and Xvariables
x = xvar;
y = fac;
ii = sum(isnan(x),2)+sum(isnan(y),2);
x = x(ii==0,:);
y = y(ii==0,:);
[~,~,r_lev] = canoncorr(x,y);   % NOTE -- Matlab demeans the series 

x = xvar_res;
y = fac_res;
ii = sum(isnan(x),2)+sum(isnan(y),2);
x = x(ii==0,:);
y = y(ii==0,:);
[~,~,r_res] = canoncorr(x,y);   % NOTE -- Matlab demeans the series 

fprintf(fileID,'  Variables included in VAR: \n');
for i = 1:size(splot.name,2);
   a1 = char(splot.name(i));
   fprintf(fileID,['     ' a1 '\n']);  
end;
fprintf(fileID,'  Sample Period: %4i:Q%1i',est_par.smpl_par.nfirst);
fprintf(fileID,'-%4i:Q%1i\n',est_par.smpl_par.nlast);
fprintf(fileID,'  Number of Factors: %3i \n',est_par.fac_par.nfac.unobserved);
fprintf(fileID,'    Canonical Correlations for Levels of X and Factor \n      ');
prtmat_comma(r_lev,fileID,'%5.2f','\n');
fprintf(fileID,'    Canonical Correlations for VAR Residuals of X and Factor \n      ');
prtmat_comma(r_res,fileID,'%5.2f','\n');


% --- Panel C ... Use Stepwise procedure to maximize correlation of VAR residuals with Factor Residuals
% Step 1: Construct dataset with variables available over full-sample period
[istart,iend] = smpl_HO(est_par.smpl_par);
ifirst = 1;
for i = 1:size(datain.bpdata,2);
    if sum(isnan(datain.bpdata(istart:iend,i)))==0;
        if ifirst == 1;
            ifirst = 0;
            bdata = datain.bpdata(:,i);
            bnamevec = datain.bpnamevec(i);
        else;
            bdata = [bdata datain.bpdata(:,i)];
            bnamevec = [bnamevec datain.bpnamevec(i)];
        end;
    end;
end;
nseries = size(bdata,2);

% -- Stepwise procedure to find variables;
for i = 1:est_par.fac_par.nfac.unobserved;
  rvec = zeros(nseries,1);
  for j = 1:nseries;
      if i == 1;
          xvar = bdata(:,j);
      else;
          xvar = [xvar_last bdata(:,j)];
      end;
      varout_xvar = varest(xvar,est_par.var_par,est_par.smpl_par);
      xvar_res = varout_xvar.resid;
      x = xvar_res;
      y = fac_res;
      ii = sum(isnan(x),2)+sum(isnan(y),2);
      x = x(ii==0,:);
      y = y(ii==0,:);
      [~,~,r_res] = canoncorr(x,y);   
      rvec(j) = r_res(end);
  end;
  [rmax,jmax]=max(rvec);
  if i == 1;
      xvar_last = bdata(:,jmax);
      xnamevec = bnamevec(jmax);
  else; 
      xvar_last = [xvar_last bdata(:,jmax)];
      xnamevec = [xnamevec bnamevec(jmax)];
  end;
  jsave = [(0:jmax-1) (jmax+1:nseries+1)];
  jsave = jsave(2:end-1);
  bdata = bdata(:,jsave);
  bnamevec = bnamevec(jsave);
  nseries = size(bdata,2);
  [i nseries rmax]
end;

fprintf(fileID,'\n\n\n\n');
% Variables to Include in VAR
splot.name = xnamevec;
% Get Series
for i = 1:size(splot.name,2);
  a1 = char(splot.name(i));
  j = colnumber(a1,datain.bpnamevec);
  y = datain.bpdata(:,j);
  if i == 1;
      xvar = y;
  else;
      xvar = [xvar y];
  end;
end;
% Compute VAR
varout_xvar = varest(xvar,est_par.var_par,est_par.smpl_par);
xvar_res = varout_xvar.resid;
% Levels of Facots and Xvariables
x = xvar;
y = fac;
ii = sum(isnan(x),2)+sum(isnan(y),2);
x = x(ii==0,:);
y = y(ii==0,:);
[~,~,r_lev] = canoncorr(x,y);   % NOTE -- Matlab demeans the series 

x = xvar_res;
y = fac_res;
ii = sum(isnan(x),2)+sum(isnan(y),2);
x = x(ii==0,:);
y = y(ii==0,:);
[~,~,r_res] = canoncorr(x,y);   % NOTE -- Matlab demeans the series 

fprintf(fileID,'  Variables included in VAR: \n');
for i = 1:size(splot.name,2);
   a1 = char(splot.name(i));
   fprintf(fileID,['     ' a1 '\n']);  
end;
fprintf(fileID,'  Sample Period: %4i:Q%1i',est_par.smpl_par.nfirst);
fprintf(fileID,'-%4i:Q%1i\n',est_par.smpl_par.nlast);
fprintf(fileID,'  Number of Factors: %3i \n',est_par.fac_par.nfac.unobserved);
fprintf(fileID,'    Canonical Correlations for Levels of X and Factor \n      ');
prtmat_comma(r_lev,fileID,'%5.2f','\n');
fprintf(fileID,'    Canonical Correlations for VAR Residuals of X and Factor \n      ');
prtmat_comma(r_res,fileID,'%5.2f','\n');
            
    





