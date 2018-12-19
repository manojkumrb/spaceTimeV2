% BLL_ReadData - upload and transforms US quarterly Data
%
% [Data, codeData] = BLL_ReadData(trans,out,sel)
%   trans = 0 - No transformation except logarithms
%   trans = 1 - Light transformation
%   trans = 2 - Heavy transformation
%   trans = 3 - Heavy Transformation but variables in Level
%   out   = 1 - Remove outliers
%   sel   = 1 - You are uploading the Real US DB
%   sel   = 2 - Standard Large US DB
%   sel   = 3 - Large US DB  with per capita variables
%

% Original version from Mario Forni, 
% modified by Matteo Luciani (matteo.luciani@.ulb.ac.be)

function [Data, Label, Name, Dates, cd] = BLL_ReadData(trans,out,sel)

load BLL_US_DB;                                                             % Upload the Database
DB(end,:)=[]; Dates(end,:)=[];                                              % accounts for the jagged-edged structure of the data

%%% Determines which Database to Upload %%%
if sel==1;     
    Select(2,:)=[]; disp('Uploading Real US DB')
elseif sel==2; 
    Select(1,:)=[]; disp('Uploading Standard Large US DB')
elseif sel==3; 
    Select(1,:)=[]; disp('Uploading Large US DB with per capita variables')
    J=[37:48 50:57 72:5:82]; 
    j1=length(J);
    DB(:,J)=DB(:,J)./repmat(DB(:,96),1,j1);
    for ii=1:j1;
        Name{J(ii)}= [Name{J(ii)} ' Per-Capita'];
    end
    
elseif sel==4;     
    Select(1,:)=[];
    population=DB(:,97)./(DB(:,95)/100);
    A=1000000000*DB(:,[37 39 54 ])./repmat(population*1000,1,3);
    J=min(find(year(Dates)==2005));
    pop100=population/population(J);
    B=DB(:,82)./pop100;           
    Trans=[Trans Trans(:,[37 39 54 82])];
    DB=[DB A B];
    Select=[Select ones(1,4)];
    Name{133}= 'Per-Capita GDP';
    Name{134}= 'Per-Capita Investment'; 
    Name{135}= 'Per-Capita Consumption';
    Name{136}= 'Per-Capita Hours';
    label{133}= 'Yh';
    label{134}= 'Ih'; 
    label{135}= 'Ch';
    label{136}= 'Hh';
end                     


jj=find(Select==1);                                                         % Selected Variables 
x=DB(:,jj); Trans=Trans(:,jj); Label=label(jj,:); Name=Name(jj,:);          % ------------------
[T,N]=size(x);                                                              % Useful objects



%%% Transform Variables %%%
if trans==2; 
    Trans(1,:)=[]; disp('Heavy Transformation');
elseif trans==3; trans=1;
    Trans(1,:)=[]; disp('Heavy Transformation but variables in Level'); 
    Trans(Trans==2)=1; Trans(Trans==3)=5; Trans(Trans==4)=3;
else
    Trans(2,:)=[]; 
    if trans==0; disp('No transformation - just take logs')
        Trans(Trans==2)=1; Trans(Trans>2)=5; 
    else disp('Light Transformation')
    end
end                         

data=NaN(T,N);                                                              % Preallocates
data(:,Trans==1)=x(:,Trans==1);                                             % Transform Variables
data(2:T,Trans==2)=ML_diff(x(:,Trans==2));                                  % -------------------
data(2:T,Trans==3)=ML_diff(100*log(x(:,Trans==3)));                         % -------------------
data(3:T,Trans==4)=ML_diff(ML_diff(100*log(x(:,Trans==4))));                % -------------------
data(:,Trans==5)=100*log(x(:,Trans==5));                                    % -------------------
Data=data(trans+1:T,:); Dates(1:trans)=[];
cd=zeros(N,1); cd(Trans==2|Trans==3)=1; cd(Trans==4)=2;

 
if out==1; Data = removeoutliers(Data); end

            %%% ========================= %%%
            %%% == ------------------- == %%%
            %%% == - Remove Outliers - == %%%
            %%% == - By Mario Forni  - == %%%
            %%% == ------------------- == %%%
            %%% ========================= %%%
            
function cleaneddata = removeoutliers(data)
kk=0;
for i = 1:size(data,2);
    iqr_data = iqr(data(:,i));                      % Interquantile Range
    amd_tdata = abs(data(:,i) - median(data(:,i))); % Distance from Median Value
    a  = find(amd_tdata >= 6*iqr_data);             % Identify Outliers    
    if ~isempty(a)
        for j = 1:length(a);
            kk=kk+1; out(kk,:)=[i a(j)];
            if a(j)>1; data(a(j),i) = median(data(max(1,a(j)-5):a(j)-1,i));
            else         data(a(j),i) =  median(data(2:6, i)); end
        end        
    end
end
% fprintf('Number of Detected Outliers'); disp(kk)
% disp(out);
cleaneddata = data;