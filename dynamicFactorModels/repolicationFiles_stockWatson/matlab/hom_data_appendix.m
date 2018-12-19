% hom_data_appendx.m
% Summarize some information for data appendix tables
% 

clear all;
small = 1.0e-10;
big = 1.0e+6;  

% -- File Directories  
outdir = '/Users/mwatson/Dropbox/Hom_Factor/ddisk/Matlab/out/';
figdir = '/Users/mwatson/Dropbox/Hom_Factor/ddisk/Matlab/fig/';
matdir = '/Users/mwatson/Dropbox/Hom_Factor/ddisk/Matlab/mat/'; 

% ----------- Sample Period, Calendars and so forth
[dnobs_m,calvec_m,calds_m] = calendar_make([1959 1],[2014 12],12);  % Monthly Calendar variables 
[dnobs_q,calvec_q,calds_q] = calendar_make([1959 1],[2014 4],4);    % Quarterly Caleendar variables

% -- Load Data
load_data=1;
  % Demeaning Parameters
  i_demean = 1;  % 0 Do Nothing
               % 1 Eliminate low-frequency by local Demeaning
               % 2 Eliminate low-frequency trend by full-sample demeaning;
    
  bw_bw = 100;   % Bi-Weight Parameter for local demeaning
datain_all;

% Easy Names for variables used here ;
bpinclcode = datain.bpinclcode;
bpcatcode = datain.bpcatcode;
bplabvec_short = datain.bplabvec_short;
bplabvec_long = datain.bplabvec_long;
bpdata_raw = datain.bpdata_raw;
bptcodevec = datain.bptcodevec;
bpoutliervec = datain.bpoutliervec;
calds = datain.calds;


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

outfile_name = [outdir 'Table_DataDescriptionSummary.out'];
fileID = fopen(outfile_name,'w');
fprintf(fileID,' ~Category~Number of series~Number of series used for factor estimation \n'); 
for i = 1:size(reordervec,2);
    fprintf(fileID,'(%-2i)~',i);
    tmp =char(catvec(i));
    fprintf(fileID,[tmp '~']);
    icat = reordervec(i);
    ii = floor(bpcatcode) == icat;
    nii = sum(ii);  % Total Number of Series in This Category;
    fprintf(fileID,'%-3i~',nii);
    jj = ii.*incl_ind;
    njj = sum(jj);  % Number of Series in this category used for Factor Estimation;
    fprintf(fileID,'%-3i \n',njj);
end;
fprintf(fileID,'~~~\n');
fprintf(fileID,'~Total~');
nii = size(bpinclcode,1);
fprintf(fileID,'%-3i~',nii);
jj = incl_ind;
njj = sum(jj);
fprintf(fileID,'%-3i \n',njj);

outfile_name = [outdir 'Table_DataDescriptionDetailed.out'];
fileID = fopen(outfile_name,'w');
fprintf(fileID,'~Name~Description~Smpl~Trans~Outlier~FacEst \n'); 
inumber = 0;
for i = 1:size(reordervec,2);
    tmp =char(catvec(i));
    fprintf(fileID,['~~' tmp '~~~ \n']);
    icat = reordervec(i);
    ii = floor(bpcatcode) == icat;
    tmp_bplabvec_short = bplabvec_short(1,ii==1);
    tmp_bplabvec_long = bplabvec_long(1,ii==1);
    tmp_bptcodevec=bptcodevec(ii==1);
    tmp_bpoutliervec=bpoutliervec(ii==1);
    tmp_bpinclcode = bpinclcode(ii==1);
    tmp_bpdata_raw = bpdata_raw(:,ii==1);
    for ij = 1:sum(ii);
        inumber = inumber+1;
        fprintf(fileID,'%-3i~',inumber);
        tmp = char(tmp_bplabvec_short(1,ij));
        fprintf(fileID,[tmp '~']);
        tmp = char(tmp_bplabvec_long(1,ij));
        fprintf(fileID,[tmp '~']);
        tmp = tmp_bpdata_raw(:,ij);
        tmp = calds(isnan(tmp)==0,:);
        fprintf(fileID,'%4i:Q%1i-%4i:Q%1i~',[tmp(1,:) tmp(end,:)]);
        fprintf(fileID,'%1i~',tmp_bptcodevec(ij));
        fprintf(fileID,'%1i~',tmp_bpoutliervec(ij));
        fprintf(fileID,'%1i \n',(tmp_bpinclcode(ij)==1));
    end;
end;

% Tabulate number of series available over complete sample period
ii = 0;
for i = 1:size(bpdata_raw,2);
    y = bpdata_raw(:,i);
    if sum(isnan(y)) == 0;
        ii = ii+1;
    end;
end;
fprintf('Number of series without missing values: %3i \n',ii);

