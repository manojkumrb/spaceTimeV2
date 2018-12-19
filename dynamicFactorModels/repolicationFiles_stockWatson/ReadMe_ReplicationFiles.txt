Replication files for Stock and Watson's "Factor Models and Structural Vector Autoregrssions in Macroeconomics", March 9, 2016


Data:  All data are in the Excel file hom_fac_1.xlsx

The data are described in the Data Appendix, and (somewhat) internally documented in the Excel data file.

Listed below are the Matlab files used to produce the various tables and figures in the paper.  These programs can be modified for use in other projects.  Feel free to do so. Please read the file SW_DFM_PROGRAMMING_NOTES.TXT for useful information about these programs.

Table 1:  hom_data_appendix.m    This program produces the results shown in Table 1 and also produces a version of the data description table in the Data Appendix.

Table 2, panel (a):  hom_descriptive_statistics_real_variables.m
Table 2, panel (b) and (c):  hom_descriptive_statistics_all_variables.m

Table 3, columns (a):  Tabulate_rsquared_static_factors.m
Table 3, columns (b):  Tabulate_variance_decomps.m  

Table 4: hom_stability.m

Table 5: hom_var_approx.m

Table 6: Tabulate_variance_decomps_1985_2014.m

SDFM, FAVAR, SVAR models (results reported in Table 7-8 and Figures 8-11)
  These programs must be run BEFORE the tables and figures and generated
		Oil Price Exogenous: svar_price_exogenous.m, favar_price_exogenous.m, and sdfm_price_exogenous.m
		Kilian Identification: svar_kilian.m, favar_kilian.m, and sdfm_kilian.m
	
Table 7: irf_oil_price_exogenous.m and irf_kilian.m 

Table 8: structural_shocks_correlation.m

Figure 1: plot_fitted_values.m

Figure 2: gain_plots.m

Figure 3: the values are computed in  hom_descriptive_statistics_real_variables. The bar graph is constructed in Scree_plots.xlsx

Figure 4  hom_figure_4.m

Figure 5: hom_figure_5.m

Figure 6a:  the values are computed in  hom_descriptive_statistics_all_variables. The bar graph is constructed in Scree_plots.xlsx

Figure 6b: figure_6b.m

Figure 7:  Oil_common_factor_plot.m

Figure 8: irf_oil_price_exogenous.m

Figures 9-11:  irf_kilian.m