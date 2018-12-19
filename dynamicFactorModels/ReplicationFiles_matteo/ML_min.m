% ML_min - Absolute min/max of a matrix of whathever dimension (it accounts for NaN)
% x=ML_min(x,type)

% Written by Matteo Luciani (matteo.luciani@ulb.ac.be)

function x=ML_min(x,type)
x(isinf(x))=NaN;
if nargin==1; type=1; end 

J=max(size(size(x)));                   % number of dimensions of matrix x

if type==1; for jj=1:J; x=nanmin(x); end   % take the minimum
else        for jj=1:J; x=nanmax(x); end   % take the maximum
end


