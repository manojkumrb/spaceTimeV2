% ML_GoodPolar - Draw a rotation matrix which satisfies soome sign restrictions
% 
% normrest - normalization restrictions 
%     rest - sign restrictions
%   ndraws - number of draws of the rotation matrix
%  horizon - number of horizon on which sign restrictions are imposed
%   maxrot -  maximum number of rotation matrices satisfying the restrictions
% 

% Written by Mario Forni
% Modified and commented by Matteo Luciani (matteo.luciani@ulb.ac.be)

function GR = ML_GoodPolar(B, normrest, rest, ndraws, horiz, maxrot)

if nargin==5; maxrot=ndraws; end
q = size(B,2);                                                              % numer of shocks
npar = q*(q-1)/2;                                                           % number of parameters
GR = ones(q,q,0);
for j=1:ndraws;
    angles = rand(npar,1)*2*pi;                                             % draws npar angles
    rotationmat = RotationMatrix(q, angles);                                % construct the rotation matrix
    effects0 = diag(B(abs(normrest),:,3)*rotationmat);                      % normalization restrictions on the sign at impact
    correctsigns = diag(sign(effects0).*sign(normrest));                    % normalization restrictions on the sign at impact
    rotationmat = rotationmat*correctsigns;                                 % invert signs of the rotation matrix if necessary
    for k = 1:horiz
        effects = sign(B(:,:,k)*rotationmat);                                                   % sign of the IRF at each horizon
        if diag(effects(rest{k}(:,1),abs(rest{k}(:,2))))'*sign(rest{k}(:,2)) < size(rest{k},1)  % Verify if sign restrictions are verified            
            k = k-1;
            break
        end        
    end
    if k == horiz(end); GR = cat(3,GR,rotationmat); end;                    % If restrictions are verified store the rotation matrix
    if size(GR,3)==maxrot; break; end;
end


function R = RotationMatrix(q, theta)
R = eye(q);
M = R;
k = 1;
for j = 1 : q - 1
    for i = j + 1 : q
        M(j , j) = cos(theta(k));
        M(j , i) = sin(theta(k));
        M(i , j) = - sin(theta(k));
        M(i , i) = cos(theta(k));
           R = R*M;
        M = eye(q);
        k = k + 1;
    end
end
