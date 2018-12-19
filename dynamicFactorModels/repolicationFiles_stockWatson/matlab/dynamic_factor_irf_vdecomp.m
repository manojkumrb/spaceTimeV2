function out = dynamic_factor_irf_vdecomp(fac_est_out, est_par, decomp_par, bptcodevec)

% Compute variance decomposition
% 
% Input:
%   fac_est_out: output from fac_est_ls_full.m function
%                contains least squares estimates of factor model
%   est_par: estimation parameters
%   decomp_par: variance decomposition parameters
%       hor: horizon for variance decomposition
%       varcum: 1 for cumulative variance of first i shocks
%       cancor: 1 for ordering shocks by cannonical correlations
%   bptcodevec: vector of codes denoting units of each time series
%


% EXTRACT INPUTS
est_data      = fac_est_out.est_data;
fac_est       = fac_est_out.fac;
lam_mat       = fac_est_out.lam_mat;
uar_coef_mat  = fac_est_out.uar_coef_mat;
uar_ser_mat   = fac_est_out.uar_ser_mat;
varout        = fac_est_out.varout;
n_varlag      = est_par.var_par.nlag;
n_uarlag      = est_par.n_uarlag;
[dnobs, nfac] = size(fac_est_out.fac);
n_series      = size(fac_est_out.lam_mat,1);
hor           = decomp_par.hor;
varcum        = decomp_par.varcum;
cancor        = decomp_par.cancor;
ntmin         = est_par.lambda.nt_min;    % minimum number of obs for regression

% If Cancor == 1, transform G so that dynamic factors are ordered to coincide with canonical correlation of Factor VAR shocks and Residuals from regression of observables onto lags of factor;
if cancor == 1;  
  % Transform G so shocks are ordered as dynamic shocks 
  % Construct lags of factors and residuals for est_data
  x = [ones(dnobs,1), lagmatrix(fac_est,1:n_varlag)];
  est_data_res = NaN(size(est_data));
  trend = (1:1:size(est_data,1))';
  for is = 1:size(est_data,2);
    tmp = [trend est_data(:,is) x];
    tmp = packr(tmp);
    if size(tmp,1) >= ntmin;
      ii = tmp(:,1);
      y = tmp(:,2);
      z = tmp(:,3:end);
      b = (z'*z)\(z'*y);
      e = y - z*b;
      est_data_res(ii,is) = e;   % Residual from regressions of data on lags of factors
    end;
  end;
   
  % Find data series available for full-sample period
  [istart,iend] = smpl_HO(est_par.smpl_par);
  istart_res = istart + n_uarlag;
  fac_r = varout.resid(istart_res:iend,:);
  y_tmp = est_data_res(istart_res:iend,:);
  y_r = [];
  for is = 1:size(y_tmp,2);
    tmp = y_tmp(:,is);
    if sum(isnan(tmp)) == 0;
        y_r = [y_r tmp];
    end;
  end;
  
  % Compute canonical variables and construct dynamic factor shocks
  % accordingly
  [~,B,~] = canoncorr(y_r,fac_r);
   % Scale G1 by DF from VAR 
    DF_m = (size(fac_r,1) - 1);
    DF_var = (size(fac_r,1) - nfac*n_varlag - 1);
    G1 = sqrt(DF_m/DF_var)*inv(B');
    varout.coef.G(1:nfac,1:nfac) = G1;

end;

% Construct IRFs for 1-SD cholesky shocks 
coef_y_fac = varout.coef; 
coef_y_fac.Q = lam_mat*coef_y_fac.Q;
imp_y_fac_mat = varirf_ss(coef_y_fac,0,hor);        % irf to factor shocks

% Compute IRF wrt to own shock
imp_y_u_mat = NaN(n_series,hor);
coef_y_u.Q = [1, zeros(1,n_uarlag-1)];
for is = 1:n_series
    coef_y_u.M = [uar_coef_mat(is,:); eye(n_uarlag-1),zeros(n_uarlag-1,1)];
    coef_y_u.G = [uar_ser_mat(is); zeros(n_uarlag-1,1)];
    imp_y_u_mat(is,:) = varirf_ss(coef_y_u,0,hor);   % irf to own shocks
end;

% Transform IRFs so they are units of levels of series
for is = 1:n_series;
    tc = bptcodevec(is);
    tmp = units_to_levels(imp_y_u_mat(is,:)', tc);
    imp_y_u_mat(is,:) = tmp';
    for i_shock = 1:nfac;
        tmp = units_to_levels(reshape(imp_y_fac_mat(is,i_shock,:),hor,1), tc);
        imp_y_fac_mat(is,i_shock,:) = tmp';
    end;
end;

% Scale Impulse Responses by standard deviations
imp_y_fac_mat_scl = NaN*zeros(n_series,nfac,hor);
for i = 1:nfac;
    imp_y_fac_mat_scl(:,i,:) = imp_y_fac_mat(:,i,:)/varout.coef.G(i,i);
end;

% Compute Variance Components for each Shock
tmp = imp_y_u_mat.^2;
vcomp_y_u_mat = cumsum(tmp,2);
vcomp_y_fac_mat = imp_y_fac_mat.^2;
if varcum == 1;
    vtotal_y_fac_mat = NaN(n_series,nfac,hor); 
else;
    vtotal_y_fac_mat = NaN(is,hor);
end;
for is = 1:n_series;
    tmp = reshape(vcomp_y_fac_mat(is,:,:),nfac,hor);
    tmp = cumsum(tmp,2);
    vcomp_y_fac_mat(is,:,:) = tmp;
    if varcum == 1
        vtotal_y_fac_mat(is,:,:) = cumsum(tmp);
    else;
        vtotal_y_fac_mat(is,:) = sum(tmp);
    end;
end;
if varcum == 1
    vfrac_y_fac_mat = NaN(is,nfac,hor);
    for is = 1:n_series
        for ii = 1:nfac
            tot = vtotal_y_fac_mat(is,ii,:);
            tot = reshape(tot,1,hor);
            frac = tot./(tot + vcomp_y_u_mat(is,:));
            vfrac_y_fac_mat(is,ii,:) = frac;
        end;
    end;
else;
    vtotal_y_mat = vtotal_y_fac_mat + vcomp_y_u_mat;   % total variance by forecast horizon
    vfrac_y_fac_mat = vtotal_y_fac_mat./vtotal_y_mat;
    vfrac_y_comp_mat = NaN(n_series,nfac,hor);
    for i = 1:nfac;
        tmp = reshape(vcomp_y_fac_mat(:,i,:),n_series,hor);
        tmp = tmp./vtotal_y_mat;
        vfrac_y_comp_mat(:,i,:) = tmp;
    end;
end;

% Compute Structural Shocks using VAR residuals and G matrix
% G*eps = resid
G = varout.coef.G(1:nfac,1:nfac);
% Impose unit coefficients on diagonal of G
% GS*eps_s = resid
SG = diag(diag(G));
GS = G*inv(SG);
% eps_s = inv(GS)*resid
GSI = inv(GS);
eps_reducedform = varout.resid;
eps_structural = eps_reducedform*GSI';

% SAVE OUTPUT
out.eps_structural = eps_structural;        % Structural Errors Given G matrix 
out.imp_y_fac_mat     = imp_y_fac_mat;      % irf to factor shocks
out.imp_y_fac_mat_scl = imp_y_fac_mat_scl;  % irf to factor shocks, scaled by standard deviation
out.vcomp_y_fac_mat   = vcomp_y_fac_mat;    % variance component from factors
out.vcomp_y_u_mat     = vcomp_y_u_mat;      % variance component from own shocks
out.vtotal_y_fac_mat  = vtotal_y_fac_mat;   % total variance from factors
out.vfrac_y_fac_mat   = vfrac_y_fac_mat;    % fraction from factor shocks
if varcum ~= 1
    out.vfrac_y_comp_mat = vfrac_y_comp_mat;   % fraction 
end;

end