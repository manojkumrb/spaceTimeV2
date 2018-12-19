% Construct Desriptive Statistics for "real" variable data set
% 11-22-2015, mww


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
datain_real;      % datain_all reads in the full dataset .. all variables, etc. saved in datain.xx

% --- Parameters used for Esimation 
  % Calendar, sample, parameters
  est_par.smpl_par.nfirst = [1959 3];         % start date
  est_par.smpl_par.nlast  = [2014 4];         % end date
  est_par.smpl_par.calvec = datain.calvec;    % calendar 
  est_par.smpl_par.nper   = 4;                % number of periods per year

  % Factor analysis parameters
  est_par.fac_par.nt_min                  = 20;     % min number of obs for any series used to est factors
  est_par.fac_par.tol                     = 10^-8;  % precision of factor estimation (scaled by by n*t)

% VAR parameters for factors
est_par.var_par.nlag   = 4;    % number of lags

% Factor Parameters
nfac_max = 11;

% Extract Estimation Sample 
est_data = datain.bpdata(:,datain.bpinclcode==1);
est_namevec = datain.bpnamevec(:,datain.bpinclcode==1);


% --- Full Sample Results ---
% Compute Various Statistics for estimating the number of static and dynamic factors
nfac_out = est_nfac(est_data,nfac_max,est_par);

% Summarize and save results
% Static factors
kvec = (1:nfac_max)';
trace_r2 = ones(nfac_max,1)- nfac_out.st.ssr/nfac_out.st.tss;
marg_r2 = NaN*ones(nfac_max,1);
marg_r2(1) = trace_r2(1);
marg_r2(2:end) = trace_r2(2:end)-trace_r2(1:end-1);
ah_er = NaN(nfac_max,1);
ah_er(1:nfac_max-1) = marg_r2(1:nfac_max-1)./marg_r2(2:nfac_max);
fprintf(fileID,'Descriptive Statistics for determining the number of factors \n\n');
fprintf(fileID,'  Sample Period: %4i:Q%1i',est_par.smpl_par.nfirst);
fprintf(fileID,'-%4i:Q%1i\n',est_par.smpl_par.nlast);
fprintf(fileID,'  Static factor statistics\n');
fprintf(fileID,'    Nobs = %8i \n',nfac_out.st.nobs);
fprintf(fileID,'    Nbar = %5.2f \n',(nfac_out.st.nobs/nfac_out.st.nt));
fprintf(fileID,'K, trace R2, marginal r2, BN-ICP, AH-ER \n');
fprintf(fileID,'0,0,0,0\n');
for i = 1:nfac_max;
    fprintf(fileID,'%2i, ', kvec(i));
    fprintf(fileID,'%5.3f, ', trace_r2(i));
    fprintf(fileID,'%5.3f, ', marg_r2(i));
    fprintf(fileID,'%5.3f, ', nfac_out.st.bn(i));
    fprintf(fileID,'%5.3f \n', ah_er(i));
 end;
 
 % AW criterion for dynamic factors
 fprintf(fileID,'\n  AW criteria for dynamic factors\n');
 fprintf(fileID,'      Rows: Dynamic factors \n');
 fprintf(fileID,'      Cols: Static factors  \n');
 fprintf(fileID,',');
 prtmat_comma(kvec',fileID,'%2i','\n');
 for i = 1:nfac_max;
     fprintf(fileID,'%2i,',i);
     prtmat_comma(nfac_out.dy.aw(i,:),fileID,'%5.3g','\n');
 end;
     

% ------- Repeat for First half of sample period
est_par.smpl_par.nfirst = [1959 3];         % start date
est_par.smpl_par.nlast  = [1983 4];         % end date
 
% Compute Various Statistics for estimating the number of static and dynamic factors
nfac_out = est_nfac(est_data,nfac_max,est_par);

% Summarize and save results
% Static factors
kvec = (1:nfac_max)';
trace_r2 = ones(nfac_max,1)- nfac_out.st.ssr/nfac_out.st.tss;
marg_r2 = NaN*ones(nfac_max,1);
marg_r2(1) = trace_r2(1);
marg_r2(2:end) = trace_r2(2:end)-trace_r2(1:end-1);
ah_er = NaN(nfac_max,1);
ah_er(1:nfac_max-1) = marg_r2(1:nfac_max-1)./marg_r2(2:nfac_max);
fprintf(fileID,'Descriptive Statistics for determining the number of factors \n\n');
fprintf(fileID,'  Sample Period: %4i:Q%1i',est_par.smpl_par.nfirst);
fprintf(fileID,'-%4i:Q%1i\n',est_par.smpl_par.nlast);
fprintf(fileID,'  Static factor statistics\n');
fprintf(fileID,'    Nobs = %8i \n',nfac_out.st.nobs);
fprintf(fileID,'    Nbar = %5.2f \n',(nfac_out.st.nobs/nfac_out.st.nt));
fprintf(fileID,'K, trace R2, marginal r2, BN-ICP, AH-ER \n');
fprintf(fileID,'0,0,0,0\n');
for i = 1:nfac_max;
    fprintf(fileID,'%2i, ', kvec(i));
    fprintf(fileID,'%5.3f, ', trace_r2(i));
    fprintf(fileID,'%5.3f, ', marg_r2(i));
    fprintf(fileID,'%5.3f, ', nfac_out.st.bn(i));
    fprintf(fileID,'%5.3f \n', ah_er(i));
 end;
 
 % AW criterion for dynamic factors
 fprintf(fileID,'\n  AW criteria for dynamic factors\n');
 fprintf(fileID,'      Rows: Dynamic factors \n');
 fprintf(fileID,'      Cols: Static factors  \n');
 fprintf(fileID,',');
 prtmat_comma(kvec',fileID,'%2i','\n');
 for i = 1:nfac_max;
     fprintf(fileID,'%2i,',i);
     prtmat_comma(nfac_out.dy.aw(i,:),fileID,'%5.3g','\n');
 end;
  
 
% Repeat for second half of sample period
est_par.smpl_par.nfirst = [1984 1];         % start date
est_par.smpl_par.nlast  = [2014 4];         % end date
 
% Compute Various Statistics for estimating the number of static and dynamic factors
nfac_out = est_nfac(est_data,nfac_max,est_par);

% Summarize and save results
% Static factors
kvec = (1:nfac_max)';
trace_r2 = ones(nfac_max,1)- nfac_out.st.ssr/nfac_out.st.tss;
marg_r2 = NaN*ones(nfac_max,1);
marg_r2(1) = trace_r2(1);
marg_r2(2:end) = trace_r2(2:end)-trace_r2(1:end-1);
ah_er = NaN(nfac_max,1);
ah_er(1:nfac_max-1) = marg_r2(1:nfac_max-1)./marg_r2(2:nfac_max);
fprintf(fileID,'Descriptive Statistics for determining the number of factors \n\n');
fprintf(fileID,'  Sample Period: %4i:Q%1i',est_par.smpl_par.nfirst);
fprintf(fileID,'-%4i:Q%1i\n',est_par.smpl_par.nlast);
fprintf(fileID,'  Static factor statistics\n');
fprintf(fileID,'    Nobs = %8i \n',nfac_out.st.nobs);
fprintf(fileID,'    Nbar = %5.2f \n',(nfac_out.st.nobs/nfac_out.st.nt));
fprintf(fileID,'K, trace R2, marginal r2, BN-ICP, AH-ER \n');
fprintf(fileID,'0,0,0,0\n');
for i = 1:nfac_max;
    fprintf(fileID,'%2i, ', kvec(i));
    fprintf(fileID,'%5.3f, ', trace_r2(i));
    fprintf(fileID,'%5.3f, ', marg_r2(i));
    fprintf(fileID,'%5.3f, ', nfac_out.st.bn(i));
    fprintf(fileID,'%5.3f \n', ah_er(i));
 end;
 
 % AW criterion for dynamic factors
 fprintf(fileID,'\n  AW criteria for dynamic factors\n');
 fprintf(fileID,'      Rows: Dynamic factors \n');
 fprintf(fileID,'      Cols: Static factors  \n');
 fprintf(fileID,',');
 prtmat_comma(kvec',fileID,'%2i','\n');
 for i = 1:nfac_max;
     fprintf(fileID,'%2i,',i);
     prtmat_comma(nfac_out.dy.aw(i,:),fileID,'%5.3g','\n');
 end;
  
