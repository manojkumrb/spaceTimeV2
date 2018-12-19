% ML_Rotation - Build a Rotation Matrix R s.t. 1) R*R' = eye = R'*R; 2) det(R) = 1
%
% [R] = ML_Rotation(theta,q);
%   R     - is q x q, it is the product of q(q-1)/2 2-dimensional rotation matrices
%   theta - q angles
%   q     - number of Shock to rotate (optional) it genrates random q angles
%  

% Written by Matteo Luciani (matteo.luciani@ulb.ac.be)

function R = ML_Rotation(theta,q)

if nargin==2; theta=2*pi*rand((q*(q-1)/2),1); end;                          % Pick initial values
    
r=size(theta,1);
q=(1+sqrt(1+8*r))/2;
kk=1;
R=eye(q);

for ii=1:q-1;
    for jj=ii+1:q;
        temp=eye(q);
        temp(ii,ii)=cos(theta(kk));
        temp(ii,jj)=sin(theta(kk));
        temp(jj,ii)=-sin(theta(kk));
        temp(jj,jj)=cos(theta(kk));          
        R=R*temp(:,:);
        kk=kk+1;
    end;
end;
