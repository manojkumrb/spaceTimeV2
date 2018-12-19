% Construct Results for Figure 6b


% --- Preliminaries --- 
clear all;
small = 1.0e-10;
big = 1.0e+6;  

% -- File Directories  
outdir = '/Users/mwatson/Dropbox/Hom_Factor/ddisk/Matlab/out/';
figdir = '/Users/mwatson/Dropbox/Hom_Factor/ddisk/Matlab/fig/';
matdir = '/Users/mwatson/Dropbox/Hom_Factor/ddisk/Matlab/mat/';

% --- File Name for output 
fname = [outdir 'balanced_panel.xls'];

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

% Calendar, sample, parameters
est_par.smpl_par.nfirst = [1959 3];         % start date
est_par.smpl_par.nlast  = [2014 4];         % end date
est_par.smpl_par.calvec = datain.calvec;    % calendar
est_par.smpl_par.nper   = 4;                % number of periods a year
    
% Extract Estimation Sample 
est_data = datain.bpdata(:,datain.bpinclcode==1);
est_namevec = datain.bpnamevec(:,datain.bpinclcode==1);

% Preliminaries
smpl_par           = est_par.smpl_par;

% Find indices variables without missing values
ns_all = size(est_data,2);
[istart, iend] = smpl_HO(smpl_par);
tmp = (1:1:ns_all)';
tmp = [tmp est_data(istart:iend,:)'];
tmp = packr(tmp);
bp_incl = tmp(:,1)';

% Carry out Calculations
% ----- Full Sample ---
smpl_par.nfirst = [1959 3];         % start date
smpl_par.nlast  = [2014 4];         % end date


% Estimate factors with unbalanced panel by LS -- standardize data first
  % Sample period
  [istart, iend] = smpl_HO(smpl_par);
  istart = max(istart,1);
  iend = min(iend,size(est_data,1));
  
  % Estimate Factors 
  xdata = est_data(istart:iend,bp_incl);
  nt = size(xdata,1);
  ns = size(xdata,2);
   
  % Mean and Standard Deviation
  xmean = mean(xdata)';   % mean 
  xstd = std(xdata)';     % std 
  xbal = (xdata - repmat(xmean',nt,1))./repmat(xstd',nt,1);   % standardized data
  
  % Compute Total Sum of Squares
  tss = 0;
  nobs = 0;
  for is = 1:ns;
      tmp = xbal(:,is);     % select series
      tmp = tmp(isnan(tmp)==0);  % drop NaN
      tss = tss+sum(tmp.^2);     % add to tss
      nobs = nobs+size(tmp,1);   % add to n*T
  end;
  
  % Estimate factors using balanced panel
  [coef,score,latent]=princomp(xbal);
  fa = score;
  
  ssr_mat = NaN(ns,ns);
  for nfac = 1:ns;
      x = fa(:,1:nfac);
      for is = 1:ns;
          y = xbal(:,is);
          lam = x\y;
          yhat = x*lam;
          e = y - yhat;
          ssr_mat(is,nfac) = sum(e.^2);
      end;
  end;
  
  tr_ssr = sum(ssr_mat,1);
  r2_vec_fs = (tss-tr_ssr')/tss;
  
  
  % ----- First Half ---
  % Estimate factors with unbalanced panel by LS -- standardize data first
  % Sample period
  smpl_par.nfirst = [1959 3];         % start date
  smpl_par.nlast  = [1984 4];         % end date

  [istart, iend] = smpl_HO(smpl_par);
  istart = max(istart,1);
  iend = min(iend,size(est_data,1));
  
  % Estimate Factors 
  xdata = est_data(istart:iend,bp_incl);
  nt = size(xdata,1);
  ns = size(xdata,2);
   
  % Mean and Standard Deviation
  xmean = mean(xdata)';   % mean 
  xstd = std(xdata)';     % std 
  xbal = (xdata - repmat(xmean',nt,1))./repmat(xstd',nt,1);   % standardized data
  
  % Compute Total Sum of Squares
  tss = 0;
  nobs = 0;
  for is = 1:ns;
      tmp = xbal(:,is);     % select series
      tmp = tmp(isnan(tmp)==0);  % drop NaN
      tss = tss+sum(tmp.^2);     % add to tss
      nobs = nobs+size(tmp,1);   % add to n*T
  end;
  
  % Estimate factors using balanced panel
  [coef,score,latent]=princomp(xbal);
  fa = score;
  
  ssr_mat = NaN(ns,ns);
  for nfac = 1:ns;
      x = fa(:,1:nfac);
      for is = 1:ns;
          y = xbal(:,is);
          lam = x\y;
          yhat = x*lam;
          e = y - yhat;
          ssr_mat(is,nfac) = sum(e.^2);
      end;
  end;
  
  tr_ssr = sum(ssr_mat,1);
  r2_vec_p1 = (tss-tr_ssr')/tss;
  
   % ----- Second Half ---
  % Estimate factors with unbalanced panel by LS -- standardize data first
  % Sample period
  smpl_par.nfirst = [1985 1];         % start date
  smpl_par.nlast  = [2014 4];         % end date

  [istart, iend] = smpl_HO(smpl_par);
  istart = max(istart,1);
  iend = min(iend,size(est_data,1));
  
  % Estimate Factors 
  xdata = est_data(istart:iend,bp_incl);
  nt = size(xdata,1);
  ns = size(xdata,2);
   
  % Mean and Standard Deviation
  xmean = mean(xdata)';   % mean 
  xstd = std(xdata)';     % std 
  xbal = (xdata - repmat(xmean',nt,1))./repmat(xstd',nt,1);   % standardized data
  
  % Compute Total Sum of Squares
  tss = 0;
  nobs = 0;
  for is = 1:ns;
      tmp = xbal(:,is);     % select series
      tmp = tmp(isnan(tmp)==0);  % drop NaN
      tss = tss+sum(tmp.^2);     % add to tss
      nobs = nobs+size(tmp,1);   % add to n*T
  end;
  
  % Estimate factors using balanced panel
  [coef,score,latent]=princomp(xbal);
  fa = score;
  
  ssr_mat = NaN(ns,ns);
  for nfac = 1:ns;
      x = fa(:,1:nfac);
      for is = 1:ns;
          y = xbal(:,is);
          lam = x\y;
          yhat = x*lam;
          e = y - yhat;
          ssr_mat(is,nfac) = sum(e.^2);
      end;
  end;
  
  tr_ssr = sum(ssr_mat,1);
  r2_vec_p2 = (tss-tr_ssr')/tss;
  
  % Plot
  trend = (1:1:ns)';
  plot(trend,r2_vec_fs,'- k','LineWidth',3);
  hold on;
    plot(trend,r2_vec_p1,'-- r','LineWidth',3);
    plot(trend,r2_vec_p2,': b','LineWidth',3);
  hold off;
  xlim([0 60]);
  ylim([0 1.0]);
  legend('1959-2014','1959-1985','1985-2014');
  legend('Location','SouthEast');
  ax = gca;
  ax.FontSize = 20;
  figname = [figdir 'Figure_6b'];
  savefig(figname);
  
  