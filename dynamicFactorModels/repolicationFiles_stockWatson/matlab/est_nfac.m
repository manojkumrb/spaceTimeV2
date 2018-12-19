function nfac_out = est_nfac(est_data,nfac_max,est_par)
% 11/21/2015, Paul Ho and Mark Watson
% 
% -- INPUT --
% est_data:  Data used for estimating factors  (nt x ns)
% nfac_max:  Max number of factors
% est_par: estimation parameters
%     A. smpl_par: sampling parameters
%         1. calvec: Calendar vector
%         2. nper:   number of periods per year
%         3. nfirst: First observation for estimation
%         4. nlast:  Last observation for estimation
%     B. fac_par: factor estimation parameters
%         1. lambda_constraints: These are linear constraints on the lambda parameters
%                                Each row of lambda_constraints corresponds to a single constraint
%                                .. first column: "i" ... which row of lambda is constrained;
%                                .. call this row lam(i)
%                                The constraint is then R*lam(i)' = r
%                                R is given in cols 2-nfac+1 of lambda_constraints
%                                r is given in last column of lambda_constraints
%                                set = 1 if not used;
%         2. nt_min:  Minimum number of observations for any series used to estimate factors
%         3. tol: Precision of estimate
%                 Loop terminates when objective function changes by less than tol/(nt*ns)
%         4. nvar_lag = Number of VAR Lags for Factors in dynamic factor calculations;
%   
% -- OUTPUT --
%  A. st: Static factor estimation
%      1. bn  (nfac_max X 1) : Bai-Ng criterion
%      2. ssr (nfac_max X 1) : sum of square residuals
%      3. r2  (ns * nfac_max): r-squared
%      4. tss (1 X 1): Total Sum of Squares
%      5. nobs(1 X 1): Total Number of Obs using for estimation
%      6. nt (1 X 1):  Total number of time periods used in estimation
%  B. dy: Dynamic factor estimation
%      1. bn  (nfac_max X nfac_max)     : Amengual-Watson criterion
%                                           row - number of dynamic factors
%                                           col - number of static factors
%      2. ssr (nfac_max X nfac_max)     : sum of square residuals
%                                           row - number of dynamic factors
%                                           col - number of static factors
%      3. r2  (ns X nfac_max X nfac_max): r-squared
%                                           dim 1 - series
%                                           dim 2 - number of dynamic factors
%                                           dim 3 - number of static factors


% ESTIMATION
bn     = NaN(nfac_max,1);                  aw     = NaN(nfac_max);
ssr_st = NaN(nfac_max,1);                  ssr_dy = NaN(nfac_max,nfac_max);
r2_st  = NaN(size(est_data,2), nfac_max);  r2_dy  = NaN(size(est_data,2), nfac_max, nfac_max);
est_par_tmp = est_par;
est_par_tmp.fac_par.lambda_constraints_est  = 1;   % no constraints on lambda 
for nfac = 1:nfac_max;
    % static factors
    est_par_tmp.fac_par.nfac.unobserved = nfac;
    est_par_tmp.fac_par.nfac.observed = 0;
    est_par_tmp.fac_par.nfac.total = nfac;
    lsout = factor_estimation_ls(est_data, est_par_tmp);
    bn(nfac)      = bai_ng(lsout);
    ssr_st(nfac)  = lsout.ssr;
    r2_st(:,nfac) = lsout.r2vec;
    % dynamic factors
    awout = amengual_watson(lsout,est_data,est_par_tmp)';
    aw(1:nfac,nfac)      = awout.aw;
    ssr_dy(1:nfac,nfac)  = awout.ssr;
    r2_dy(:,1:nfac,nfac) = awout.r2;
end;

% SAVE OUTPUT
nfac_out.st.bn  = bn;         
nfac_out.st.ssr = ssr_st;     
nfac_out.st.r2  = r2_st;      
nfac_out.st.tss = lsout.tss;
nfac_out.st.nobs = lsout.nobs;
nfac_out.st.nt = lsout.nt;
nfac_out.dy.aw  = aw;
nfac_out.dy.ssr = ssr_dy;
nfac_out.dy.r2  = r2_dy;

end