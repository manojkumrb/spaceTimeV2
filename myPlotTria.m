% plot tria mesh
tri=load('tria.txt');
tri=tri+ones(size(tri));
coOrdinate= load('half.xyz');
trimesh(tri,coOrdinate(:,1),coOrdinate(:,2),coOrdinate(:,3));