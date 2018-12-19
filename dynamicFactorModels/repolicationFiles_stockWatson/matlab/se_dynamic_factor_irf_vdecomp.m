function out = se_dynamic_factor_irf_vdecomp(datain,fac_est_out,est_par,decomp_par,n_rep)

% Compute standard errors for irfs and variance decomps using (Gaussian) parametric bootstrap
% 
% Input:
%   fac_est_out: output from fac_est_ls_full_HO.m function
%                contains least squares estimates of factor model
%   est_par: estimation parameters
%   decomp_par: variance decomposition parameters
%       hor: horizon for variance decomposition
%       varcum: 1 for cumulative variance of first i shocks
%       cancor: 1 for ordering shocks by cannonical correlations
%   bptcodevec: vector of codes denoting units of each time series
%   n_rep = number of bootstap simulations

% Construct Gaussian parametric bootstrap draws ;

% Bootstrap Esimation Parameters
% ... Factors (including observed factors are simulated
bs_est_par = est_par;

bs_est_par.fac_par.w = 1;

% ------ Easy names for paramters ----
[dnobs, n_fac] = size(fac_est_out.fac);
n_series      = size(datain.bpdata,2);
uar_coef_mat = fac_est_out.uar_coef_mat;
uar_ser_mat = fac_est_out.uar_ser_mat;
M = fac_est_out.varout.coef.M;
G = fac_est_out.varout.coef.G;
n_uarlag = est_par.n_uarlag;
n_varlag = est_par.var_par.nlag;
n_initial = 100;      % number of discarded initial values for VAR and AR simulations;

% Matrices Used below;
bpdata_NaN = isnan(datain.bpdata);
Mu = zeros(n_uarlag,n_uarlag);
if n_uarlag > 1;
  Mu(2:end,1:end-1) = eye(n_uarlag-1);
end;
Gu = zeros(n_uarlag,1);

% --- Matrices for saving sum and sum^2 of draws
sum_btstrp_imp_y_fac_mat_scl = zeros(n_series,n_fac,decomp_par.hor);
sum_squared_btstrp_imp_y_fac_mat_scl = zeros(n_series,n_fac,decomp_par.hor);
sum_btstrp_vfrac_y_comp_mat = zeros(n_series,n_fac,decomp_par.hor);
sum_squared_btstrp_vfrac_y_comp_mat = zeros(n_series,n_fac,decomp_par.hor);

for irep = 1:n_rep;
  t1 = cputime;
  % ---- Generate Bootstap Data
  % Generate Factor;
  f_btstrp = zeros(n_fac*n_varlag,n_initial+dnobs);      
  for t = 2:dnobs+n_initial;        % bootstrap sample of {f(t)}
    f_btstrp(:,t) =  M*f_btstrp(:,t-1) + G*randn(size(G,2),1);
  end
  f_btstrp = f_btstrp(1:n_fac,n_initial+1:end);  % Note this is nfac X dnobs;
  % Generate u
  u_btstrp = NaN(dnobs,n_series);
  for i = 1:n_series;
    Mu(1,:) = uar_coef_mat(i,:);
    Gu(1,1) = uar_ser_mat(i);
    u = zeros(n_uarlag,n_initial+dnobs);
    for t = 2:dnobs+n_initial;
      u(:,t) =  Mu*u(:,t-1) + Gu*randn(1,1);
    end;
    u_btstrp(:,i) = u(1,n_initial+1:end)';
  end;
  % Generate Y
  y_btstrp = (fac_est_out.lam_mat*f_btstrp)' + u_btstrp;
  % Replace values with missing values as in original dataset
  y_btstrp(bpdata_NaN==1) = NaN;
  % Save oberved factors
  if bs_est_par.fac_par.nfac.observed > 0;
     bs_est_par.fac_par.w = f_btstrp(1:bs_est_par.fac_par.nfac.observed,:)';
  end;
  
  % Compute VD and IRF
  btstrp_fac_est_out = factor_estimation_ls_full(y_btstrp,datain.bpinclcode,bs_est_par);                             % estimation
  btstrp_irf_vdecomp_out = dynamic_factor_irf_vdecomp(btstrp_fac_est_out,bs_est_par,decomp_par,datain.bptcodevec);   % irf variance decomposition
  sum_btstrp_imp_y_fac_mat_scl = sum_btstrp_imp_y_fac_mat_scl + btstrp_irf_vdecomp_out.imp_y_fac_mat_scl;
  sum_squared_btstrp_imp_y_fac_mat_scl = sum_squared_btstrp_imp_y_fac_mat_scl+ btstrp_irf_vdecomp_out.imp_y_fac_mat_scl.^2;
  sum_btstrp_vfrac_y_comp_mat = sum_btstrp_vfrac_y_comp_mat + btstrp_irf_vdecomp_out.vfrac_y_comp_mat;
  sum_squared_btstrp_vfrac_y_comp_mat =  sum_squared_btstrp_vfrac_y_comp_mat + btstrp_irf_vdecomp_out.vfrac_y_comp_mat.^2; 
  
  if irep == 5;
    t2 = cputime-t1;
    fprintf('Carry out simulations to compute standard errors \n');
    fprintf('  Number of simulations: %6i \n',n_rep);
    fprintf('  Seconds per simulation: %6.3f \n',t2);
    fprintf('  Total time (in seconds) required %10.2f \n',t2*n_rep);
  end;
end;

mean_btstrp_imp_y_fac_mat_scl = sum_btstrp_imp_y_fac_mat_scl/n_rep;
var_btstrp_imp_y_fac_mat_scl = sum_squared_btstrp_imp_y_fac_mat_scl/n_rep - mean_btstrp_imp_y_fac_mat_scl.^2;
mean_btstrp_vfrac_y_comp_mat = sum_btstrp_vfrac_y_comp_mat/n_rep;
var_btstrp_vfrac_y_comp_mat =  sum_squared_btstrp_vfrac_y_comp_mat/n_rep - mean_btstrp_vfrac_y_comp_mat.^2;
se_btstrp_imp_y_fac_mat_scl = sqrt(var_btstrp_imp_y_fac_mat_scl); 
se_btstrp_vfrac_y_comp_mat = sqrt(var_btstrp_vfrac_y_comp_mat);

out.mean_imp_y_fac_mat_scl = mean_btstrp_imp_y_fac_mat_scl;
out.mean_vfrac_y_comp_mat = mean_btstrp_vfrac_y_comp_mat;
out.se_imp_y_fac_mat_scl = se_btstrp_imp_y_fac_mat_scl;
out.se_vfrac_y_comp_mat = se_btstrp_vfrac_y_comp_mat;

end