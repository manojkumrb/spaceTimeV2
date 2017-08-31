%% to interpolate eigen vectors
% as fem structure is used for element indentification.. we create a new
% mesh from the triangilation developed in CGAL for the key nodes
% import eigen vectors corresponding to the selected key nodes
% load all the nodes corresponding to dense mesh in 'points'

close all; clear ;
dbstop if error;

run 'C:\Users\babu_m\Documents\GitHub\source\initFiles.m';

fem2=femInit();
fileName{1}=strcat(inpFiles,'\hinge_25_90.inp');%\top_hat.inp');%\halo.inp');%\hinge.inp');%
% fem=importMesh(fem,fileName{1});%'top_hat.inp'); % use only for top hat
fem2=importMultiMesh(fem2, fileName);
fem2=femPreProcessing(fem2);
totalNodes=size(fem2.xMesh.Node.Coordinate,1);
domainID=1;


points=load('fullHingeMeshCoordi.mat'); % dense mesh
points=points.hingeCoordi;              %
V=load('eigVec25_90.mat');%eigVec72nodes.mat');% % corresponding eigen vector
V=V.V;

nodedom=fem2.Domain(domainID).Node;
keyNodes=fem2.xMesh.Node.Coordinate(nodedom,:);
interpEigVec=zeros(size(points,1),size(V,2));

denseMeshId = dsearchn(points,keyNodes);

vGlobal		= zeros(size(points,1),size(V,2));

% % setting eigen vectors for global nodeIDs
for i= 1:size(V,1)
	vGlobal(denseMeshId(i),:)=V(i,:);
end


Pp=zeros(size(points,1),3) ;
elem=zeros(size(points,1),1) ;					% element insde which the node is present
baryCoeff=zeros(size(points,1),3) ;
errCount=0;
errNode=zeros(size(points,1),1) ;
eigVec=zeros(size(points,1),size(V,2)) ;


for i=1:size(points,1)
	
	% STEP 1: GLOBAL SEARCH
	temp=keyNodes;
	P0=points(i,:);
	temp(:,1)=temp(:,1)-P0(1);
	temp(:,2)=temp(:,2)-P0(2);
	temp(:,3)=temp(:,3)-P0(3);
	
	dik=sqrt(sum(temp.^2,2));
	[dm, mid]=min(dik);
	midl=mid;						% keeping local index as well
	mid=nodedom(mid);				% getting global node ID
	%----------------------
	% STEP 2: LOCAL SEARCH
	try
		flag=false;
		
		[flag,~,~,Pp(i,:),elem(i),baryCoeff(i,:)]=getProjectionElement(fem2,mid,P0,50);
		
	catch
		if flag==false;
			% 			warning('node outside triangilation: setting eigen vector equal to the closest node')
			Pp(i,:)=keyNodes(mid,:) ;
			
			errCount=errCount+1;
			errNode(errCount)=i;
		else
			warning('unkonowm error check!')
			sprintf('at node %i',i)
		end
	end
	
end

errNode(errCount+1:end)=[]; % removing empty entries

meshplotAxisDefined(fem2,1);
hold all
plot3(points(errNode,1),points(errNode,2),points(errNode,3),'*')
title('unProjected Nodes','fontweight','bold','fontsize',12);
ax=gca;
ax.Clipping='off';
hold off



%% setting to nearest node
for i=1:size(points,1)
	
	
	if elem(i)==0   % if node is not inside the triangulation then consider the nearest node
		interpEigVec(i,:)=V(mid,:);		
	else              % interpolation of eigen values according to wikle and cressie
		
		tria=fem2.xMesh.Element(elem(i)).Element;
		baryCoord=baryCoeff(i,:);
		interpEigVec(i,:)=baryCoord(1).*V(tria(1),:)+baryCoord(2).*V(tria(2),:)+baryCoord(3).*V(tria(3),:);
		
	end
	
end

%% insted of using nearest node use nearest within the mesh
withInMesh=find(elem);
coordWithin=points(withInMesh,:);

for i=1:size(points,1)
	if elem(i)==0   % if node is not inside the triangulation then consider the nearest node
		
		localID=dsearchn(coordWithin,points(i,:));
		golbaID=withInMesh(localID);
		interpEigVec(i,:)=interpEigVec(golbaID,:);
		
	end
end


% %% checking points and corresponding elements
meshplotAxisDefined(fem2,domainID)
hold all;
ax=gca;
ax.Clipping='off';
for i=1:20
	
	golbaID=withInMesh(i);
	tria=fem2.xMesh.Element(elem(golbaID)).Element;
	plot3(points(golbaID,1),points(golbaID,2),points(golbaID,3),'*')
	plot3(keyNodes(tria,1),keyNodes(tria,2),keyNodes(tria,3),'o')
	
end


%% testing other interpolation
% nInterp=3000;
% flag=zeros(nInterp,1);%zeros(size(points,1),1);
% interpEigVec2=zeros(nInterp,size(V,2));%zeros(size(points,1),size(V,2));
% for i=1:size(V,2)
% 	fem2.Post.Interp.InterpVariable='user';
% 	fem2.Post.Interp.Pm=Pp(1:nInterp,:);	
% 	fem2.Sol.UserExp = V(:,i);
% 
% 	[~, interpEigVec2(:,i), flag(:)]=getInterpolationData_fast(fem2);
% 
% end

% contourDomainPlot(fem2,1,V(:,1),0)
% hold all
% plot3(points(withInMesh(1:20),1),points(withInMesh(1:20),2),points(withInMesh(1:20),3),'*')

