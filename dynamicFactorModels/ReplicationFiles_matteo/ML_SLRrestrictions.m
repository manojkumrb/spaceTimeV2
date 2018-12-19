% ML_SLRrestrictions - Short and Long Run zero restrictions
% 
% H=ML_SLRrestrictions(NN,theta,B1,K0,K1)
% 
%   NN = variables on which restrictions are imposed
%   theta = Impulse response functions (N x q x s)
%   B1 = long run IRF
%   K0 = shor run restrictions
%   K1 = Long Run restrictions
% 

% Written by Matteo Luciani (matteo.luciani@ulb.ac.be)

function [H,eflag,x]=ML_SLRrestrictions(NN,theta,B1,K0,K1)

q=size(theta,2);                                                            % number of shocks
C0=theta(NN,:,1); C1=B1(NN,:);                                              % Selecting multipliers  
x0=2*pi*rand((q*(q-1)/2),1);                                                % Pick initial values
options=optimset('Display','off','MaxFunEvals',1000);                       % options for maximization
[x,~,eflag]=fsolve(@(x) ML_identification(x,C0,C1,K0,K1),x0,options);       % Estimate the angles for the Rotation Matrix
% [x,~,eflag]=lsqnonlin(@(x) ML_identification(x,C0,C1,K0,K1),...
%     x0,zeros(q,1),2*pi*ones(q,1),options);% Estimate the angles for the Rotation Matrix

H=ML_Rotation(x);                                                           % Rotation Matrix
v=diag(C0*H); ii=find(sign(v)<0); H(:,ii)=-H(:,ii);                         % Normalization


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F=ML_identification(x,C0,C1,K0,K1)
r=length(x);
F=zeros(r,1);
H = ML_Rotation(x);                                                         % Rotation Matrix 
B0=C0*H; B1=C1*H;                                                           % Structural Multiplier 
k0=size(K0,1);  if K0(1,1)==0; k0=0; end                                    % Number of Short-run Restrictions
for ii=1:k0;   F(ii)=B0(K0(ii,1),K0(ii,2)); end;                            % Short run restrictions
for ii=k0+1:r; F(ii)=B1(K1(ii-k0,1),K1(ii-k0,2)); end;                      % Long-Run Restrictions


function F=ML_identification2(x,C0,C1,K0,K1)
r=length(x);
F=zeros(r,1);
H = ML_Rotation(x);                                                         % Rotation Matrix 
B0=C0*H; B1=C1*H;                                                           % Structural Multiplier 
k0=size(K0,1);  if K0(1,1)==0; k0=0; end                                    % Number of Short-run Restrictions
for ii=1:k0;   F(ii)=B0(K0(ii,1),K0(ii,2)); end;                            % Short run restrictions
for ii=k0+1:r; F(ii)=B1(K1(ii-k0,1),K1(ii-k0,2)); end;                      % Long-Run Restrictions

F=sum(F);