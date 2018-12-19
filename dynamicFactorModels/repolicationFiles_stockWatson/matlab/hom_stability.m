% Construct Desriptive Statistics for real variable data set
% 11-22-2015, mww


% --- Preliminaries --- 
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
datain_all;        % datain_1 reads in the full dataset

% Calendar, sample, parameters
est_par.smpl_par.calvec = datain.calvec;    % calendar
est_par.smpl_par.nper   = 4;                % number of periods a year

% Factor analysis parameters
%nfac.unobserved = 4;
nfac.unobserved = 8;
nfac.observed = 0;
nfac.total = nfac.unobserved+nfac.observed;
est_par.fac_par.nfac = nfac;
est_par.fac_par.lambda_constraints_est  = 1;      % no constraints on lambda
est_par.fac_par.nt_min                  = 20;     % min number of obs for any series used to est factors
est_par.lambda.nt_min                   = 40;     % min number of obs for any series used to estimate lamba, irfs, etc.
est_par.fac_par.tol                     = 10^-8;  % precision of factor estimation (scaled by by n*t)
est_par.fac_par.nvar_lag                = 4;      % Number of VAR lags for factors in computign dynamic factor model        

% Extract Estimation Sample 
est_data = datain.bpdata(:,datain.bpinclcode==1);
est_namevec = datain.bpnamevec(:,datain.bpinclcode==1);

% Full-sample factors
est_par.smpl_par.nfirst = [1959 3];         % start date
est_par.smpl_par.nlast  = [2014 4];         % end date
lsout_fs = factor_estimation_ls(est_data, est_par);

% First half-sample factors
est_par.smpl_par.nfirst = [1959 3];         % start datei
est_par.smpl_par.nlast  = [1984 4];         % end date
lsout_s1 = factor_estimation_ls(est_data, est_par);

% Second half-sample factors
est_par.smpl_par.nfirst = [1985 1];         % start date
est_par.smpl_par.nlast  = [2014 4];         % end date
lsout_s2 = factor_estimation_ls(est_data, est_par);

% Sample Period Indicators
ismpl_fs = smpl(datain.calvec,[1959 3],[2014 4],4);
ismpl_s1 = smpl(datain.calvec,[1959 3],[1984 4],4);
ismpl_s2 = smpl(datain.calvec,[1985 1],[2014 4],4);

% Data series
bpdata = datain.bpdata;
ns = size(bpdata,2);
dnobs = size(bpdata,1);
chow_vec = NaN(ns,1);
qlr_vec = NaN(ns,1);
r_s1 = NaN(ns,1);
r_s2 = NaN(ns,1);

% Loop through series and compute statistics;
for i = 1:ns;
    y = bpdata(:,i);
    y1 = y(ismpl_s1==1);
    y2 = y(ismpl_s2==1);
    % Check and carry out analysis if there are >= 80 quarters in each half
    % sample
    n1 = sum(isnan(y1)==0);
    n2 = sum(isnan(y2)==0);
    nmin = min([n1;n2]);
    if nmin >= 80;
        x = [lsout_fs.fac];
        y = bpdata(:,i);
        x_2 = x.*repmat(ismpl_s2,1,size(x,2));
        z = [x x_2];
        tmp = [y z];
        tmp = packr(tmp);
        yp = tmp(:,1);
        zp = tmp(:,2:end); 
        nma = 6;
        ikern = 1;
        [betahat,vbeta,se_beta,ser,rbarsq] = hac(yp,zp,nma,ikern);
        i1 = size(x,2);
        b1 = betahat(i1+1:end);
        v1 = vbeta(i1+1:end,i1+1:end);
        chow = b1'*(inv(v1))*b1;
        chow_vec(i) = chow;
        
        % Compute QLR Statistic
        x = [lsout_fs.fac];
        y = bpdata(:,i);
        z = x;
        tmp = [y z];
        tmp = packr(tmp);
        yp = tmp(:,1);
        zp = tmp(:,2:end); 
        nma = 6;
        [lm, lmr, lsbreak] = qlra(yp,1,zp,0.15,nma);
        qlr_vec(i) = lmr; 
        
        % Compute Fitted values
        yfit_fs = NaN(dnobs,1);
        yfit_s1 = NaN(dnobs,1);
        yfit_s2 = NaN(dnobs,1);
        trend = (1:1:dnobs)';
        
        % Full Sample
        x = [lsout_fs.fac];
        y = bpdata(:,i);
        tmp = [trend y x];
        tmp = packr(tmp);
        ii_fs = tmp(:,1);
        yp = tmp(:,2);
        xp = tmp(:,3:end);
        b = xp\yp;
        yp_fit = xp*b;
        yfit_fs(ii_fs) = yp_fit;
        
        % First sub-sample
        x = [lsout_s1.fac];
        y = bpdata(:,i);
        tmp = [trend y x];
        tmp = packr(tmp);
        ii_s1 = tmp(:,1);
        yp = tmp(:,2);
        xp = tmp(:,3:end);
        b = xp\yp;
        yp_fit = xp*b;
        yfit_s1(ii_s1) = yp_fit;
        
        % Second sub-sample
        x = [lsout_s2.fac];
        y = bpdata(:,i);
        tmp = [trend y x];
        tmp = packr(tmp);
        ii_s2 = tmp(:,1);
        yp = tmp(:,2);
        xp = tmp(:,3:end);
        b = xp\yp;
        yp_fit = xp*b;
        yfit_s2(ii_s2) = yp_fit;
        
        % Compute Correlations
        tmp = packr([yfit_fs yfit_s1]);
        y1 = tmp(:,1);
        y2 = tmp(:,2);
        y1 = y1 - mean(y1);
        y2 = y2 - mean(y2);
        r_s1(i) = y1'*y2/(sqrt(y1'*y1)*sqrt(y2'*y2));
        
        tmp = packr([yfit_fs yfit_s2]);
        y1 = tmp(:,1);
        y2 = tmp(:,2);
        y1 = y1 - mean(y1);
        y2 = y2 - mean(y2);
        r_s2(i) = y1'*y2/(sqrt(y1'*y1)*sqrt(y2'*y2));
    end;
end;

chow_df = nfac.total;
chow_cv = chi2inv([0.99 0.95 0.90],chow_df);
chow_vecp = chow_vec(isnan(chow_vec)==0);
chow_rej = repmat(chow_vecp,1,3) > repmat(chow_cv,size(chow_vecp,1),1);

qlr_vecp = qlr_vec(isnan(qlr_vec)==0);
if nfac.total == 8;
    qlr_cv = [3.57 2.98 2.69]*8;   % CVs from SW text
    qlr_rej = repmat(qlr_vecp,1,3) > repmat(qlr_cv,size(qlr_vecp,1),1);
end;
if nfac.total == 4;
    qlr_cv = [5.12 4.09 3.59]*4;   % CVs from SW text
    qlr_rej = repmat(qlr_vecp,1,3) > repmat(qlr_cv,size(qlr_vecp,1),1);
end;

r_s1p = packr(r_s1);
r_s2p = packr(r_s2);
pct = [0.05 0.25 0.50 0.75 0.95];
r_s1_percentile = pctile(r_s1p,pct);
r_s2_percentile = pctile(r_s2p,pct);

% Save Results

% --- File Name for output 
outfile_name = [outdir 'Stability_Analysis_' num2str(nfac.unobserved) '.out'];
fileID = fopen(outfile_name,'w');

fprintf(fileID,'Results from Stability Tests \n');
fprintf(fileID,'Number of Factors: %2i \n',nfac.total);
fprintf(fileID,'Number of Series with 80 or more quarters in both samples: %3i \n',size(chow_vecp,1));
fprintf(fileID,'Chow Break tests for lambda (break in 1984Q4)  \n');
tmp = mean(chow_rej);
fprintf(fileID,'  Fraction rejected at 1, 5 and 10 percent levels: %5.2f  %5.2f  %5.2f \n',tmp);
fprintf(fileID,'QLR Break tests for lambda (break in 1984Q4)  \n');
tmp = mean(qlr_rej);
fprintf(fileID,'  Fraction rejected at 1, 5 and 10 percent levels: %5.2f  %5.2f  %5.2f \n',tmp);
fprintf(fileID,' \n\n Percentiles of distribution of correlations between full-sample and split sample fitted values \n');
prtmat_comma(pct,fileID,'%5.2f','\n');
fprintf(fileID,'Full sample and first sub-sample:');
prtmat_comma(r_s1_percentile',fileID,'%5.2f','\n');
fprintf(fileID,'Full sample and second sub-sample:');
prtmat_comma(r_s2_percentile',fileID,'%5.2f','\n');


% Tabulate Results by category
incl_ind = bpinclcode == 1;
% List Series in Each Category
cat1 = 'NIPA';
cat2 = 'Industrial Production';
cat3 = 'Employment and Unemployment';
cat5 = 'Orders, Inventories, and Sales';
cat4 = 'Housing Starts and Permits';
cat6 = 'Prices';
cat7 = 'Productivity and Earnings';
cat8 = 'Interest Rates';
cat9 = 'Money and Credit';
cat12 = 'International Variables';
cat10 = 'Asset Prices, Wealth, and Household Balance Sheets';
cat14 = 'Other';
cat20 = 'Oil Market Variables';

% reorder
reordervec = [1 2 3 5 4 6 7 8 9 12 10 14 20];
catvec = {cat1 cat2 cat3 cat5 cat4 cat6 cat7 cat8 cat9 cat12 cat10 cat14 cat20};
fprintf(fileID,'\n\n\n');
fprintf(fileID,'Chow Break at 5% level by category \n');
chow_cv = chi2inv(0.95,chow_df);


fprintf(fileID,' ~Category~Number of series~Number of series avaiable~rejection frac~Median_R2_1~Median_R2_2 \n'); 
ichow = isnan(chow_vec)==0;
for i = 1:size(reordervec,2);
    fprintf(fileID,'(%-2i)~',i);
    tmp =char(catvec(i));
    fprintf(fileID,[tmp '~']);
    icat = reordervec(i);
    ii = floor(bpcatcode) == icat;
    nii = sum(ii);  % Total Number of Series in This Category;
    fprintf(fileID,'%-3i~',nii);
    jj = (ii.*ichow);
    njj = sum(jj);
    fprintf(fileID,'%-3i~',njj);
    chow_vecp = chow_vec(jj==1);
    chow_rej = mean(chow_vecp > chow_cv);
    fprintf(fileID,'%5.2f~',chow_rej);
    r_s1p = r_s1(jj==1);
    r_s1_median = pctile(r_s1p,0.50);
    r_s2p = r_s2(jj==1);
    r_s2_median = pctile(r_s2p,0.50);   
    fprintf(fileID,'%5.2f~',r_s1_median);
    fprintf(fileID,'%5.2f \n',r_s2_median);
end;






