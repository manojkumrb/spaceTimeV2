% Tabulate Correlation of Structural Shocks, Various identifications
%
 
clear all;
small = 1.0e-10;
big = 1.0e+6;
rng(63287545);

% -- File Directories  
outdir = '/Users/mwatson/Dropbox/Hom_Factor/ddisk/Matlab/out/';
figdir = '/Users/mwatson/Dropbox/Hom_Factor/ddisk/Matlab/fig/';
matdir = '/Users/mwatson/Dropbox/Hom_Factor/ddisk/Matlab/mat/'; 

outfile_name = [outdir 'struc_shocks_correlation.out'];
fileID = fopen(outfile_name,'w');

% -- Read in Structural errors 
str_ouput = 'sdfm_kilian';
str_tmp = [matdir 'rslt_' str_ouput]; load(str_tmp,'rslt');
eps_sdfm_kilian = rslt.irf_vdecomp_out.eps_structural;

str_ouput = 'favar_kilian';
str_tmp = [matdir 'rslt_' str_ouput]; load(str_tmp,'rslt');
eps_favar_kilian = rslt.irf_vdecomp_out.eps_structural;

str_ouput = 'svar_kilian';
str_tmp = [matdir 'rslt_' str_ouput]; load(str_tmp,'rslt');
eps_svar_kilian = rslt.irf_vdecomp_out.eps_structural;

str_ouput = 'sdfm_oil_price_exogenous';
str_tmp = [matdir 'rslt_' str_ouput]; load(str_tmp,'rslt');
eps_sdfm_ope = rslt.irf_vdecomp_out.eps_structural;

str_ouput = 'favar_oil_price_exogenous';
str_tmp = [matdir 'rslt_' str_ouput]; load(str_tmp,'rslt');
eps_favar_ope = rslt.irf_vdecomp_out.eps_structural;

str_ouput = 'svar_oil_price_exogenous';
str_tmp = [matdir 'rslt_' str_ouput]; load(str_tmp,'rslt');
eps_svar_ope = rslt.irf_vdecomp_out.eps_structural;

% Calendar vector
calvec = rslt.datain.calds;

tmp = [calvec eps_sdfm_ope(:,1) eps_favar_ope(:,1) eps_svar_ope(:,1)];
tmp = [tmp eps_sdfm_kilian(:,1) eps_favar_kilian(:,1)  eps_svar_kilian(:,1)];
tmp = [tmp eps_sdfm_kilian(:,2) eps_favar_kilian(:,2)  eps_svar_kilian(:,2)];
tmp = [tmp eps_sdfm_kilian(:,3) eps_favar_kilian(:,3)  eps_svar_kilian(:,3)];
tmp = packr(tmp);
e = tmp(:,3:end);
d1 = tmp(1,1:2);
d2 = tmp(end,1:2);

% Standardize
nt = size(e,1);
emean = mean(e)';   % mean 
estd = std(e)';     % std 
ebal = (e - repmat(emean',nt,1))./repmat(estd',nt,1);   % standardized data
cor = ebal'*ebal/(nt-1); 


fprintf(fileID,'Correlation of Structural Shocks \n');
fprintf(fileID,'Order is \n');
fprintf(fileID,'sdfm_oil_price_exogenous \n');
fprintf(fileID,'favar_oil_price_exogenous \n');
fprintf(fileID,'svar_oil_price_exogenous \n');
fprintf(fileID,'sdfm_kilian_supply \n');
fprintf(fileID,'favar_kilian_supply \n');
fprintf(fileID,'svar_kilian_supply \n');
fprintf(fileID,'sdfm_kilian_global_demand \n');
fprintf(fileID,'favar_kilian_global_demand \n');
fprintf(fileID,'svar_kilian_global_demand \n');
fprintf(fileID,'sdfm_kilian_oil_specific_demand \n');
fprintf(fileID,'favar_kilian_oil_specific_demand \n');
fprintf(fileID,'svar_kilian_oil_specific_demand \n\n');

for i = 1:size(cor,1);
   prtmat_comma(cor(i,:),fileID,'%4.2f','\n');
end;

