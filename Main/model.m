%%                        BIMSLOPE                       %%
%%   An open source code for slope stability analysis    %%
%%      of block-in-matrix (BIM) geomaterial             %%
%%                                                       %%
%%                 Author: Qing-xiang Meng               %%
%%                Supervisor: D.V. Griffith              %%
%%                                                       %%
%%                                                       %%
%%                     Hohai University                  %%
%%                 Version 1.0 ; July 2019               %%
%%                                                       %%
clear;clc;
rand('state',0)
w1=20;s1=20;w2=20;
h1=10;h2=5;
% generate a slope model with six points
point=zeros(4,2);
point(1,:)=[0,0];
point(2,:)=[20,0];
point(3,:)=[20,20];
point(4,:)=[0,20];
visual_poly(point);
domain=polyshape(point(:,1),point(:,2));
% Generate points for particle location
dt=0.25;
n1=(20)/dt;
n2=(20)/dt;
x=linspace(0,20,n1+1);
y=linspace(0,20,n2+1);
[X,Y]=meshgrid(x,y);
n=(n1+1)*(n2+1);
pos=[reshape(X,n,1),reshape(Y,n,1)];
[in,on] = inpolygon(pos(:,1),pos(:,2),point(:,1),point(:,2));
plot(pos(in,1),pos(in,2),'.b'); hold on
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
istop=0;
for I=1:length(rgrade)-1
    if istop==1
        break;
    end
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
                if particlenum>30
                    istop=1;
                    break;
                end
                Polygon{particlenum}=polyshape(poly(:,1),poly(:,2));
                pos_all=update_pos_2d(poly,pos_all);
                asum=asum+polyarea(poly(:,1),poly(:,2));
                Tsum=Tsum+polyarea(poly(:,1),poly(:,2));
                tmpnum=0;
%                 disp(['Generate particle ',num2str(particlenum),' fraction is ',num2str(Tsum/sarea)]);
                visual_poly(poly);
            end
        else
            tmpnum=tmpnum+1;
        end
    end
    disp(['grade ',num2str(I),' is ok']);
end
axis equal
axis off