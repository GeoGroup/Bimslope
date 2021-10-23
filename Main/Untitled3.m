g_coord=vert';
g_num=[tria(:,1),tria(:,3),tria(:,2)]';
dmax=max(sqrt(disps(:,1).^2+disps(:,2).^2));
disp(['max displacemnet is ',num2str(dmax)]);
fid=fopen('slope2.vtk','wt');
fprintf(fid,'# vtk DataFile Version 3.0\n');
fprintf(fid,'Bimslope\n');
fprintf(fid,'ASCII\n');
fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
fprintf(fid,'\n');
nn=size(vert,1);
nels=length(tnum);
fprintf(fid,'POINTS  %d   float\n',nn);
for i=1:nn
    fprintf(fid,'%f %f %f \n',g_coord(1,i),g_coord(2,i),0);
end
fprintf(fid,'\n');
fprintf(fid,'CELLS %d %d\n',nels,4*nels);
for i=1:nels
    fprintf(fid,'%d %d %d %d \n',3,g_num(1,i)-1,g_num(2,i)-1,g_num(3,i)-1);
end
fprintf(fid,'\n');
fprintf(fid,'CELL_TYPES %d\n',nels);
for i=1:nels
    fprintf(fid,' 5\n');
end
fprintf(fid,'\n');
fprintf(fid,'POINT_DATA %d\n',nn);
fprintf(fid,'VECTORS disp float\n');
for i=1:nn
    fprintf(fid,'%f %f %f\n',(disps(i,1)),(disps(i,2)),0);
end
fprintf(fid,'CELL_DATA %d\n',nels);
fprintf(fid,'SCALARS group double 1\n');
fprintf(fid,'LOOKUP_TABLE table1\n');
for i =1:nels
    fprintf(fid,'%d\n',tnum(i));
end
fclose(fid);