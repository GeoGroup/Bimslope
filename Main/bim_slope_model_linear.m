%%                        BIMSLOPE                       %%
%%   An open source code for slope stability analysis    %%
%%      of block-in-matrix (BIM) geomaterial             %%
%%                                                       %%
%%                 Author: Qing-xiang Meng               %%
%%                Supervisor: Weiya Xu                   %%
%%                            D.V. Griffith              %%
%%                                                       %%
%%                     Hohai University                  %%
%%                 Version 1.0 ; July 2019               %%
%%                                                       %%
clear;clc;
rng(1)
w1=20;s1=20;w2=20;
h1=10;h2=5;
% generate a slope model with six points
point=zeros(6,2);
point(1,:)=[0,0];
point(2,:)=[w1+s1+w2,0];
point(3,:)=[w1+s1+w2,h2];
point(4,:)=[w1+s1,h2];
point(5,:)=[w1,h1+h2];
point(6,:)=[0,h1+h2];
% visual_poly(point);
domain=polyshape(point(:,1),point(:,2));
% Generate points for particle location
dt=0.25;
n1=(w1+s1+w2)/dt;
n2=(h1+h2)/dt;
x=linspace(0,w1+s1+w2,n1+1);
y=linspace(0,h1+h2,n2+1);
[X,Y]=meshgrid(x,y);
n=(n1+1)*(n2+1);
pos=[reshape(X,n,1),reshape(Y,n,1)];
[in,on] = inpolygon(pos(:,1),pos(:,2),point(:,1),point(:,2));
% plot(pos(in,1),pos(in,2),'.b');
pos_all=pos(in,:);
sarea=polyarea(point(:,1),point(:,2));
%
addpath('C:\Users\m-q-x\OneDrive\Code\rand_comp\2d');
dr=0.25;   eta=0.0; n1=8;  n2=12;  ext=1.25;
rot=-1;distance=0.1;
rgrade=[2.5 1.5 0.75 0.5];fraction=50;frac=[65 25 10];
% fraction=30;
Polygon=cell(1000,1);
particlenum=0;
Tsum=0;
tmpnum1=0;
asum=0;atotal=0;
for I=1:length(rgrade)-1
    color=[ 1/length(rgrade)*I 0 0];
    rmax=rgrade(I);rmin=rgrade(I+1);
    atotal=atotal+sarea*fraction/100*frac(I)/100;
    
    tmpnum=0;
    particle_loop=0;
    JJ=0;
    while asum<atotal && tmpnum<20000
        JJ=JJ+1;
        tmpnum=tmpnum+1;
        p = particle_generate_polygon_2d(rmin,rmax,dr,eta,n1,n2,ext,rot);
        pos=pos_all(unidrnd(size(pos_all,1)),:);
        poly=[p(:,1)+pos(1),p(:,2)+pos(2)];
        poly2=addthick(poly,distance);
%         visual_poly(p);
        in = inpolygon(poly2(:,1),poly2(:,2),point(:,1),point(:,2));
        if sum(in)==size(poly2,1)
            collision_particle=poly_collistion(poly2,Polygon,particlenum);
            if collision_particle==1
                %break;
                tmpnum=tmpnum+1;
            else
                particlenum=particlenum+1;
                Polygon{particlenum}=polyshape(poly(:,1),poly(:,2));
                pos_all=update_pos_2d(poly2,pos_all);
                asum=asum+polyarea(poly(:,1),poly(:,2));
                Tsum=Tsum+polyarea(poly(:,1),poly(:,2));
                tmpnum=0;
                disp(['Generate particle ',num2str(particlenum),' fraction is ',num2str(Tsum/sarea)]);
%                 visual_poly(poly);
            end
        else
            tmpnum=tmpnum+1;
        end
    end
    disp(['grade ',num2str(I),' is ok']);
end
% disp(['fraction is ',num2str(Tsum/sarea)]);
% plot(pos_all(:,1),pos_all(:,2),'.');

% %---------------------------------------------- create geom and mesh
addpath('..\mesh2d');

% edge1=zeros(size(point,1),3);
node=zeros(10000,2);
edge=zeros(10000,3);
np=0;ne=0;
nod1=point;label=0;
edg1 = edge_from_vertice(point,np,label);
np=np+size(point,1);
ne=ne+size(point,1);
node(1:np,:)=nod1;
edge(1:np,:)=edg1;
for i=1:particlenum
    nod=Polygon{i}.Vertices;
    edg=edge_from_vertice(nod,np,1);
    nt=np+size(nod,1);
    node(np+1:nt,:)=nod;
    edge(np+1:nt,:)=edg;
    np=np+size(nod,1);
end
node=node(1:np,:);
edge=edge(1:np,:);



% np=0;
%     nod1=point;
%     edg1=edge_from_vertice(point,np,0);
%    nod2=Polygon{1}.Vertices;
%    np=size(point,1);
%     edg2=edge_from_vertice(nod2,np,1);
%     edg2(:,3) = +1;   
%     edg2(:,1:2) = ...
%     edg2(:,1:2);
%     edge = [edg1; edg2];
%     node = [nod1; nod2];
        
    
%-- the PART argument is a cell array that defines individu-
%-- al polygonal "parts" of the overall geometry. Each elem-
%-- ent PART{I} is a list of edge indexes, indicating which
%-- edges make up the boundary of each region.
    part{1} = [ ...
        find(edge(:,3) == 0) 
        find(edge(:,3) == 1)    
        ] ;
    part{2} = [ ...
        find(edge(:,3) == 1)
        ] ;
%     part{3} = [ ...
%         find(edge(:,3) == 2)
%         ] ;
        
    edge = edge(:,1:2) ;
    
%---------------------------------------------- do size-fun.
    hmax = +0.25;
 
   [vlfs,tlfs, ...
    hlfs] = lfshfn2(node,edge, ...
                    part) ;
    
    hlfs = min(hmax,hlfs) ;
    
   [slfs] = idxtri2(vlfs,tlfs) ;
   
%---------------------------------------------- do mesh-gen.
    hfun = @trihfn2;

   [vert,etri, ...
    tria,tnum] = refine2(node,edge,part,  [], ...
                         hfun, ...
                         vlfs,tlfs,slfs,hlfs) ;
                         
%---------------------------------------------- do mesh-opt.
   [vert,etri, ...
    tria,tnum] = smooth2(vert,etri,tria,tnum) ;
    
    figure;
    patch('faces',tria(tnum==1,1:3),'vertices',vert, ...
        'facecolor',[1.,1.,1.], ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',tria(tnum==2,1:3),'vertices',vert, ...
        'facecolor',[.9,.9,.9], ...
        'edgecolor',[.2,.2,.2]) ;
    patch('faces',tria(tnum==3,1:3),'vertices',vert, ...
        'facecolor',[.8,.8,.8], ...
        'edgecolor',[.2,.2,.2]) ;
    patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    title(['MESH-OPT.: KIND=DELFRONT, |TRIA|=', ...
        num2str(size(tria,1))]) ;

% Numerical simulation using FEM

nnode=size(vert,1);
nelement=length(tnum);
node=vert;
%% get the constraint node
% get the bottom
t=node(:,2);
s=(t<=0.1);id=1:1:nnode;
sid=s.*id';sid(sid==0)=[];
bottom=sid;
% get the left
t=node(:,1);
s=(t<=0.1);id=1:1:nnode;
sid=s.*id';sid(sid==0)=[];
left=sid;
% get the right
t=node(:,1);
s=(abs(t-w1-s1-w2)<0.1);id=1:1:nnode;
sid=s.*id';sid(sid==0)=[];
right=sid;
% form nf [node_number x_fix y_fix]
nr=ones(length(right)+length(left)+length(bottom),3);
nr(1:length(bottom),1)=bottom; nr(1:length(bottom),2:3)=0;
nr(length(bottom)+1:length(bottom)+length(left)+length(right),1)=[left;right];
nr(length(bottom)+1:length(bottom)+length(left)+length(right),2)=0;
etype=tnum;
ndim=2;nodof=2;nod=3;element='triangle' ;
nip=1;
nels=nelement;
nn=nnode;
g_coord=node';
g_num=[tria(:,1),tria(:,3),tria(:,2)]';
nf=ones(nodof,nn);
for i=1:size(nr,1)
    nf(:,nr(i,1))=nr(i,2:nodof+1);
end
nf=formnf(nf);
neq=max(nf(:));

ndof=nod*nodof; % 一个单元的自由度个数
g_g=zeros(ndof,nels); %

% ! -----------------------loop the elements to find global arrays sizes-----
for iel=1:nels
    %     if mod(iel,100000)==0
    %         disp(num2str(iel));
    %     end
    num=g_num(:,iel);
    [ g ] = num_to_g(num,nf);
    g_num(:,iel)=num;
    g_g(:,iel)=g ;
    %     kdiag = fkdiag(kdiag, g);
end
% calculate node force
gravlo=zeros(neq,1); %重力导致的节点力
[points, weights]= sample(element,nip );

KK=sparse(neq,neq);
interval=100000;
t1=mod(nels,interval);
ih=4;
I=zeros(100000000,1);J=I;K=I;idx=0;
for iel=1:nels
    if mod(iel,100000)==0
        KM=sparse(I(1:idx,:),J(1:idx,:),K(1:idx,:),neq,neq);
        KK=KK+KM;
        I=zeros(100000000,1);J=I;K=I;idx=0;
    end
        eld=zeros(ndof,1);
        dee=deemat(2000,0.25,ih);
        num=g_num(:,iel);
        coord=g_coord(:,num)';
        g=g_g(:,iel);
        km=zeros(ndof,ndof);
        % gauss integration
        for i=1:nip
            fun=shape_fun(points,i,ndim,nod);
            der=shape_der(points,i,ndim,nod);
            jac=der*coord;
            Det=det(jac);
            deriv=jac\der;
            bee=beemat( deriv,ih);
            km=km+bee'*dee*bee*Det*weights(i);
            eld(nodof:nodof:ndof)=eld(nodof:nodof:ndof)+fun(:)*Det*weights(i);
        end
         for i=1:ndof
            if g(i)~=0
                gravlo(g(i))=gravlo(g(i))-eld(i)*23.5;
                for j=1:ndof
                    if g(j)~=0
                        idx=idx+1;
                        I(idx)=g(i);J(idx)=g(j);K(idx)=km(i,j);
                    end
                end
            end
         end
end

KM=sparse(I(1:idx,:),J(1:idx,:),K(1:idx,:),neq,neq);
KK=KK+KM;
% disp('end loop');
loads=zeros(neq,1);

loads=loads+gravlo;
loads= cholmod2(KK,loads);
disp('finish mechanical calculation');
% loads(abs(loads)<1e-10)=0;
disps=zeros(nn,nodof);

for i=1:nn
    for j=1:nodof
        if nf(j,i)==0
            disps(i,j)=0;
        else
            disps(i,j)=loads(nf(j,i));
        end
    end
end

fid=fopen('slope2.vtk','wt');
fprintf(fid,'# vtk DataFile Version 3.0\n');
fprintf(fid,'Bimslope\n');
fprintf(fid,'ASCII\n');
fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
fprintf(fid,'\n');
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
    fprintf(fid,'%f %f %f\n',(disps(i,1)),abs(disps(i,2)),0);
end
fclose(fid);

% xy=node(tria(1,:),:);