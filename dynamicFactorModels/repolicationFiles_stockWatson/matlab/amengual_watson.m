function awout = amengual_watson(lsout, est_data, est_par)

% take output from fac_est_ls and compute number of dynamic factors

% PRELIMINARIES
fac = lsout.fac;             % extract estimated factors
dnobs = size(fac,1);         % number of observations
nfac_static = size(fac,2);   % number of static factors
% Sampling parameters
smpl_par = est_par.smpl_par;
nfirst   = smpl_par.nfirst;
nper     = smpl_par.nper;
% Factor estimation parameters
nvar_lag = est_par.var_par.nlag;

% Construct lags of factors and residuals for est_data
x = [ones(dnobs,1), lagmatrix(fac,1:nvar_lag)];
est_data_res = NaN*zeros(size(est_data));
trend = (1:1:size(est_data,1))';
for is = 1:size(est_data,2);
    tmp = [trend est_data(:,is) x];
    tmp = packr(tmp);
    ii = tmp(:,1);
    y = tmp(:,2);
    z = tmp(:,3:end);
    ndf = size(z,1)-size(z,2);
    if ndf >= est_par.fac_par.nt_min;   % Minimum degrees of freedom for series
      b = z\y;
      e = y - z*b;
      est_data_res(ii,is) = e;
    end;
end;

% Carry out calculations for number of dynamic factors 
nfq = nfirst(2) + nvar_lag;
ny = floor(nfq/nper);
nq = nfq-nper*ny;
est_par_res = est_par;
est_par_res.smpl_par = smpl_par;
est_par_res.smpl_par.nfirst = nfirst;
est_par_res.smpl_par.nfirst(1) = nfirst(1) + ny;
est_par_res.smpl_par.nfirst(2) = nq;
ssr = NaN(nfac_static,1);
r2  = NaN(size(est_data,2),nfac_static);
aw  = NaN(nfac_static,1);
for nfac = 1:nfac_static;
    est_par_res.fac_par.nfac.unobserved = nfac;
    est_par_res.fac_par.nfac.observed = 0;
    est_par_res.fac_par.nfac.total = nfac;
    lsout = factor_estimation_ls(est_data_res, est_par_res);
    aw(nfac)   = bai_ng(lsout);
    ssr(nfac)  = lsout.ssr;
    r2(:,nfac) = lsout.r2vec;
end


% Save output
awout.aw  = aw;
awout.ssr = ssr;
awout.r2  = r2;

end