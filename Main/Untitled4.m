fid=fopen('tslope3.vtk','wt');
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
% fprintf(fid,'POINT_DATA %d\n',nn);
fprintf(fid,'VECTORS disp2 float\n');
for i=1:nn
    fprintf(fid,'%f %f %f\n',(disps(i,1)),(disps(i,2)),0);
end
% fprintf(fid,'VECTORS stress float\n');
% for i=1:nn
%     fprintf(fid,'%f %f %f\n',nstress(1,i),nstress(2,i),nstress(3,i));
% end
% fprintf(fid,'VECTORS strain float\n');
% for i=1:nn
%     fprintf(fid,'%f %f %f\n',nstrain(1,i),nstrain(2,i),nstrain(3,i));
% end
fprintf(fid,'CELL_DATA %d\n',nels);
fprintf(fid,'SCALARS group float 1\n');
fprintf(fid,'LOOKUP_TABLE table1\n');
for i =1:nels
    fprintf(fid,'%d\n',tnum(i));
end
fprintf(fid,'SCALARS stress float 4\n');
fprintf(fid,'LOOKUP_TABLE table2\n');
for i =1:nels
    fprintf(fid,'%f %f \n',stress(1,i),stress(2,i),stress(3,i),stress(4,i));
end
fprintf(fid,'SCALARS strain float 4\n');
fprintf(fid,'LOOKUP_TABLE table3\n');
for i =1:nels
    fprintf(fid,'%f %f \n',strain(1,i),strain(2,i),strain(3,i),strain(4,i));
end
fclose(fid);