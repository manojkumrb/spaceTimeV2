% Plot IRFs from Kilian Identification SDFM, FAVAR, SVAR
%
 
clear all;
small = 1.0e-10;
big = 1.0e+6;
rng(63287545);

% -- File Directories  
outdir = '/Users/mwatson/Dropbox/Hom_Factor/ddisk/Matlab/out/';
figdir = '/Users/mwatson/Dropbox/Hom_Factor/ddisk/Matlab/fig/';
matdir = '/Users/mwatson/Dropbox/Hom_Factor/ddisk/Matlab/mat/'; 

% Fig name
fig_name_str = [figdir 'IRF_Kilian_'];

outfile_name = [outdir 'Table7_Kilian.out'];
fileID = fopen(outfile_name,'w');

% -- SDFM IRFs
str_ouput = 'sdfm_kilian';
str_tmp = [matdir 'rslt_' str_ouput]; load(str_tmp,'rslt');
datain_sdfm = rslt.datain;
irf_vdecomp_out_sdfm = rslt.irf_vdecomp_out; 
se_irf_vdecomp_out_sdfm = rslt.se_irf_vdecomp_out;
% -- favar IRFs
str_ouput = 'favar_kilian';
str_tmp = [matdir 'rslt_' str_ouput]; load(str_tmp,'rslt');
datain_favar = rslt.datain;
irf_vdecomp_out_favar = rslt.irf_vdecomp_out; 
se_irf_vdecomp_out_favar = rslt.se_irf_vdecomp_out;
% -- SDFM IRFs
str_ouput = 'svar_kilian';
str_tmp = [matdir 'rslt_' str_ouput]; load(str_tmp,'rslt');
datain_svar = rslt.datain;
irf_vdecomp_out_svar = rslt.irf_vdecomp_out; 
se_irf_vdecomp_out_svar = rslt.se_irf_vdecomp_out;


% -- Variables to Plot 
str_name = {...
'WPU0561'; ...
'OILPROD_SA'; ...
'GLOBAL_ACT'; ...
'CPIGAS'; ...
'GDPC96'; ...
'PCECC96'; ...
'FPIC96_Q'; ...
'PAYEMS'; ...
'LNS14000000'; ...
'PCECTPI'; ...
'JCXFE'; ...
'FEDFUNDS' ...
};

str_title = {...
'(a) Oil Prices'; ...
'(b) Oil Production'; ...
'(c) Global Commodity Demand Indicator'; ...
'(d) Gasoline Prices (CPI)'; ...
'(e) GDP'; ...
'(f) Consumption'; ...
'(g) Fixed Investment'; ...
'(h) Employment'; ...
'(i) Unemployment Rate'; ...
'(j) Inflation'; ...
'(k) Core Inflation'; ...
'(l) Federal Funds Rate' ...
};

str_yaxes = {...
'Percentage Points'; ...
'Percentage Points'; ...
'Standard Deviations'; ...
'Percentage Points'; ...
'Percentage Points'; ...
'Percentage Points'; ...
'Percentage Points'; ...
'Percentage Points'; ...
'Percentage Points'; ...
'Percentage Points'; ...
'Percentage Points'; ...
'Percentage Points' ...
};

scale_vec = [1 1 1 1 1 1 1 1 100 1 1 100 ]';
svar_ind = [1 1 1 0 1 0 0 1 0 1 0 1];  % = 1 if in SVAR, 0 otherwise ;
     
nimp1 = 17;
nimp2 = 13;
trend = (0:1:nimp2-1)';

nshock = 3;
for ishock = 1:nshock;
    
figure;
for i = 1:6;
    str = char(str_name(i));
    j = colnumber(str,datain_sdfm.bpnamevec);
    irf = reshape(irf_vdecomp_out_sdfm.imp_y_fac_mat_scl(j,ishock,:),1,nimp1);
    se = reshape(se_irf_vdecomp_out_sdfm.se_imp_y_fac_mat_scl(j,ishock,:),1,nimp1); 
    irf_p = irf+se;
    irf_m = irf-se;
    irf_sdfm = irf/scale_vec(i);
    irf_sdfm_p = irf_p/scale_vec(i);
    irf_sdfm_m = irf_m/scale_vec(i);
    
    irf = reshape(irf_vdecomp_out_favar.imp_y_fac_mat_scl(j,ishock,:),1,nimp1);
    se = reshape(se_irf_vdecomp_out_favar.se_imp_y_fac_mat_scl(j,ishock,:),1,nimp1); 
    irf_p = irf+se;
    irf_m = irf-se;
    irf_favar = irf/scale_vec(i);
    irf_favar_p = irf_p/scale_vec(i);
    irf_favar_m = irf_m/scale_vec(i);
    
    irf = reshape(irf_vdecomp_out_svar.imp_y_fac_mat_scl(j,ishock,:),1,nimp1);
    se = reshape(se_irf_vdecomp_out_svar.se_imp_y_fac_mat_scl(j,ishock,:),1,nimp1); 
    irf_p = irf+se;
    irf_m = irf-se;
    irf_svar = irf/scale_vec(i);
    irf_svar_p = irf_p/scale_vec(i);
    irf_svar_m = irf_m/scale_vec(i);
        
    subplot(3,2,i);
    plot(trend,irf_sdfm(1:nimp2),'- b','LineWidth',3);
    hold on;
      plot(trend,irf_sdfm_p(1:nimp2),'- b','LineWidth',1);;
      plot(trend,irf_sdfm_m(1:nimp2),'- b','LineWidth',1);
      plot(trend,irf_favar(1:nimp2),'-- r','LineWidth',3);
      if svar_ind(i) == 1;
       plot(trend,irf_svar(1:nimp2),': k','LineWidth',3);
      end;
    hold off;
    str1 = char(str_title(i));
    title(str1,'FontSize',18);
    str2 = char(str_yaxes(i));
    %ylabel(str2,'FontSize',18);
    %xlabel('Quarters','FontSize',18);
    xlim([0 12]);
    ax = gca;
    ax.FontSize = 20;  
end;
fname = [fig_name_str 'Shock_' num2str(ishock) '_PanelA'];
savefig(fname);

figure;
for i = 7:12;
    str = char(str_name(i));
    j = colnumber(str,datain_sdfm.bpnamevec);
    irf = reshape(irf_vdecomp_out_sdfm.imp_y_fac_mat_scl(j,ishock,:),1,nimp1);
    se = reshape(se_irf_vdecomp_out_sdfm.se_imp_y_fac_mat_scl(j,ishock,:),1,nimp1); 
    irf_p = irf+se;
    irf_m = irf-se;
    irf_sdfm = irf/scale_vec(i);
    irf_sdfm_p = irf_p/scale_vec(i);
    irf_sdfm_m = irf_m/scale_vec(i);
    
    irf = reshape(irf_vdecomp_out_favar.imp_y_fac_mat_scl(j,ishock,:),1,nimp1);
    se = reshape(se_irf_vdecomp_out_favar.se_imp_y_fac_mat_scl(j,ishock,:),1,nimp1); 
    irf_p = irf+se;
    irf_m = irf-se;
    irf_favar = irf/scale_vec(i);
    irf_favar_p = irf_p/scale_vec(i);
    irf_favar_m = irf_m/scale_vec(i);
    
    irf = reshape(irf_vdecomp_out_svar.imp_y_fac_mat_scl(j,ishock,:),1,nimp1);
    se = reshape(se_irf_vdecomp_out_svar.se_imp_y_fac_mat_scl(j,ishock,:),1,nimp1); 
    irf_p = irf+se;
    irf_m = irf-se;
    irf_svar = irf/scale_vec(i);
    irf_svar_p = irf_p/scale_vec(i);
    irf_svar_m = irf_m/scale_vec(i);
        
    subplot(3,2,i-6);
    plot(trend,irf_sdfm(1:nimp2),'- b','LineWidth',2);
    hold on;
      plot(trend,irf_sdfm_p(1:nimp2),': b','LineWidth',2);;
      plot(trend,irf_sdfm_m(1:nimp2),': b','LineWidth',2);
      plot(trend,irf_favar(1:nimp2),'-- r','LineWidth',2);
      if svar_ind(i) == 1;
       plot(trend,irf_svar(1:nimp2),'-. k','LineWidth',2);
      end;
    hold off;
    str1 = char(str_title(i));
    title(str1,'FontSize',18);
    str2 = char(str_yaxes(i));
    %ylabel(str2,'FontSize',18);
    %xlabel('Quarters','FontSize',18);
    xlim([0 12]);
    ax = gca;
    ax.FontSize = 20;  
end;
fname = [fig_name_str 'Shock_' num2str(ishock) '_Panelb'];
savefig(fname);

% Write VD at particular lag
% Results for Table 7

str_name2 = {...
'GDPC96'; ...
'PCECC96'; ...
'FPIC96_Q'; ...
'PAYEMS'; ...
'LNS14000000'; ...
'PCECTPI'; ...
'JCXFE'; ...
'FEDFUNDS'; ...
'MCOILBRENTEU'; ...
'OILPROD_SA'; ...
'GLOBAL_ACT'; ...
'CPIGAS' ...
};

nhor = 6;
fprintf(fileID,'Shock Number: %2i \n',ishock);
for i = 1:size(str_name2,1);
    str = char(str_name2(i));
    j = colnumber(str,datain_sdfm.bpnamevec);
    vd_favar=irf_vdecomp_out_favar.vfrac_y_comp_mat(j,ishock,nhor);
    vd_sdfm=irf_vdecomp_out_sdfm.vfrac_y_comp_mat(j,ishock,nhor);
    fprintf(fileID,[str ',']);
    fprintf(fileID,'%4.2f ,',vd_favar);
    fprintf(fileID,'%4.2f \n',vd_sdfm);
end;
fprintf(fileID,'\n\n');

end;

