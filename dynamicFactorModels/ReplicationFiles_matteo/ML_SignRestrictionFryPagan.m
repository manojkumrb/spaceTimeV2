% ML_SignRestrictionFryPagan - Implement the Fry & Pagan Correction for sign restrictions
% 
% BR=ML_SignRestrictionFryPagan(Imp,rest)
%   BR - best rotation
%  Imp - the N x q x s x ngr containing the impulse responses
% rest - a cell{k}(nv,2), by default the procedure pick the impulse at impact
% 
% Fry & Pagan "Sign Restrictions in Structural Vector Autoregressions: A
%               Critical Review", Journal of Economic Literature
% 

% Written Matteo Luciani (matteo.luciani@ulb.ac.be)

function BR=ML_SignRestrictionFryPagan(Imp,rest)

if size(Imp,4)==1; BR=1;
else   
    q=size(Imp,2);
    ngr=size(Imp,4);                                                        % Number of good rotations
    nv=rest{1}(:,1);                                                        % number of variables involved in the identifications i.e. on which are applied the restrictions
    v=[];
    for qq=1:q; temp=squeeze(Imp(nv,qq,1,:)); v = cat(1,v,temp); end        % Store the Impulse in ngr vetors of dimension nv*q by 1 
    sv=std(v,[],2);                                                         % standard deviation of impulse responses across models
    mv=median(v,2);                                                         % median of impulse responses across models
    vv=(v-mv*ones(1,ngr))./(sv*ones(1,ngr));                                % standardized impulse across models
    BR=ML_argmin(sum(vv.*vv)');                                             % identify the best rotation as in Fry Pagan
end