function lambda_constraints = lambda_construct(var_sub, var_full, R, r)

% Construct Contraints on lambda matrices

ustr = char(var_sub);
j_estdata = colnumber(ustr,var_full);
R_stack = kron(ones(length(var_sub),1),R);
r_stack = kron(ones(length(var_sub),1),r);
lambda_constraints = [kron(j_estdata,ones(size(R,1),1)),R_stack,r_stack];

end