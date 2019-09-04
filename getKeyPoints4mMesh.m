function [keyNodeCoord,denseMeshId]=getKeyPoints4mMesh(filename,domainID,points)
%
% gets the key point co-ordinates and global id from the coarse mesh
%
% init		= complete path for the initialisation file
% domainID	= domain ID of the coarse mesh to use
% filename	= file name of the coarse file if in the working folder else
% the complete path to it
% points	= the co-ordinates of the global (dense) mesh
%
% Dependencies: VRM loaded and should be in path

fem2		=femInit();
fileName{1}	=filename;%strcat(inpFiles,'\',file); % input files 
fem2		=importMultiMesh(fem2, fileName);
fem2		=femPreProcessing(fem2);
% totalNodes	=size(fem2.xMesh.Node.Coordinate,1);



nodedom=fem2.Domain(domainID).Node;
keyNodeCoord=fem2.xMesh.Node.Coordinate(nodedom,:);
denseMeshId = dsearchn(points,keyNodeCoord);

% % to check the interpolated mesh
% meshplotAxisDefined(fem2,domainID);
% hold all
% scatter3(points(denseMeshId,1),points(denseMeshId,2),points(denseMeshId,3),'*')

end