%%                        BIMSLOPE                       %%
%%   An open source code for slope stability analysis    %%
%%      of block-in-matrix (BIM) geomaterial             %%
%%                                                       %%
%%                 Author: Qing-xiang Meng               %%
%%              Supervisor: W.Y Xu & D.V. Griffith       %%
%%                                                       %%
%%                                                       %%
%%                     Hohai University                  %%
%%                 Version 1.0 ; July 2019               %%
%%                                                       %%
clear;clc;
tic;
rng('default');
rng(1);
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
%                 disp(['Generate particle ',num2str(particlenum),' fraction is ',num2str(Tsum/sarea)]);
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
    hmax = +0.5;
 
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
    

% Numerical simulation using FEM
% prop=[30,15,0,20,1.0e5,0.3;60,150,0,20,1.0e6,0.1;]';
% prop=[30,15,0,20,1.0e5,0.3;60,3000,0,20,1.0e5,0.3;]';
prop=[30,15,0,20,1.0e5,0.3
    60,300,0,30,5.0e5,0.3]';
error=1;fos_down=0;fos_up=0;iter_tot=0;
srf=1.0;istop=0;
[fail,disps,stress,strain] =fem_sim_slope(srf,vert,tria,tnum,prop);
iter_tot=iter_tot+1;
disp(['iter num ',num2str(iter_tot)]);
if fail==1
    disp('Bim slope is not safe with sf<1');
    istop=1;
else
    fos_down=1.0;
end
if istop==0
    srf=5.0;
    [fail,disps,stress,strain] =fem_sim_slope(srf,vert,tria,tnum,prop);
    iter_tot=iter_tot+1;
    disp(['iter num ',num2str(iter_tot)]);
    if fail==0
        disp('Bim slope is too safe with sf>5');
        istop=1;
    else
        fos_up=5.0;
    end
end
if istop==0
    while error>2.0e-2
        srf=(fos_down+fos_up)/2;
        [fail,disps,stress,strain] =fem_sim_slope(srf,vert,tria,tnum,prop);
        iter_tot=iter_tot+1;
        disp(['iter num ',num2str(iter_tot)]);
        if fail==0
            fos_down=srf;
        else
            fos_up=srf;
        end
        error=fos_up-fos_down;
        disp(['error is ',num2str(error)]);
    end
end
dmax=max(sqrt(disps(:,1).^2+disps(:,2).^2));
disp(['max displacemnet is ',num2str(dmax)]);
out_vtk('myslope',vert,tria,tnum,disps,stress,strain);
toc;