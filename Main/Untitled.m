% Numerical simulation using FEM
% prop=[30,15,0,20,1.0e5,0.3;60,150,0,20,1.0e6,0.1;]';
% prop=[30,15,0,20,1.0e5,0.3;60,3000,0,20,1.0e5,0.3;]';
prop=[30,15,0,20,1.0e5,0.3
    30,15,0,20,1.0e6,0.3]';
srf=1.8;
nnode=size(vert,1);
nelement=length(tnum);
node=vert;
etype=tnum;
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
% etype=tnum;
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
    if mod(iel,100000)==0 % This section can svae memory
        KM=sparse(I(1:idx,:),J(1:idx,:),K(1:idx,:),neq,neq);
        KK=KK+KM;
        I=zeros(100000000,1);J=I;K=I;idx=0;
    end
        eld=zeros(ndof,1);
        dee=deemat(prop(5,etype(iel)),prop(6,etype(iel)),ih);
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
                gravlo(g(i))=gravlo(g(i))-eld(i)*prop(4,etype(iel));
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
% loads=zeros(neq,1);
% loads=loads+gravlo;
% loads= cholmod2(KK,loads);
% disp('finish mechanical calculation');
% % loads(abs(loads)<1e-10)=0;
% disps=zeros(nn,nodof);
zero=0;two=2; dt=1.0e15;d3=3.0;d4=4.0;d6=6;one=1.0;%nprops=6;
penalty=1e20;start_dt=1.e15;nip=1;nst=4;
ndim=2;
penalty=1e20;start_dt=1.e15;nst=4;
iy=1;
np_types=2;
dt=start_dt;tol=1e-4;limit=500;
for i=1:np_types
    phi=prop(1,i);
    tnph=tan(phi*pi/180);
    phif=atan(tnph/srf(iy));
    snph=sin(phif) ;
    e=prop(5,i) ;
    v=prop(6,i);
    ddt=4*(1+v)*(1-2*v)/(e*(1-2*v+snph^2));
    if ddt<dt
        dt=ddt;
    end
end
    iters=0;
    bdylds=zeros(neq,1);
    evpt=zeros(nst,nip,nels);
    oldis=zeros(neq,1);
    iterk=0;
    %-----------------------plastic iteration loop----------------------------
    while 1

        fmax=0;
        iters=iters+1 ;
        loads=gravlo+bdylds ;
        loads= cholmod2(KK,loads);
        tensor=zeros(nst,nip,nels);
        if (iy==1 && iters==1)
            elastic=loads;
        end
        % !-----------------------check plastic convergence-------------------------
        [converged,oldis,res] = checon(loads,oldis,tol);
                iterk=iterk+1;
        disp(['iteration is  ',num2str(iterk),'  error is ',num2str(res) ]);
        if iters==1
            converged=0;
        end
        if converged==1 || iters==limit
            bdylds=zeros(neq,1);
        end
        % !-----------------------go round the Gauss Points -----------------------
        for iel=1:nels %elements_3
%             if etype(iel)==1
            bload=zeros(ndof,1);
            phi=prop(1,etype(iel));
            tnph=tan(phi*pi/180);
            phif=atan(tnph/srf(iy))*180/pi ;
            psi=prop(3,etype(iel));
            tnps=tan(psi*pi/180) ;
            psif=atan(tnps/srf(iy))*180/pi;
            cf=prop(2,etype(iel))/srf(iy) ;
            e=prop(5,etype(iel)) ;
            v=prop(6,etype(iel)) ;
            dee = deemat(e,v,nst);
            rload=zeros(ndof,1);% 
            num=g_num(:,iel);
            coord=g_coord(:,num)';
            g=g_g(:,iel);
            eld=zeros(length(g),1);
            for ii=1:length(g) %eld=loads(g) 
                if g(ii)~=0
                    eld(ii)=loads(g(ii));
                else
                    eld(ii)=0;
                end
            end
            for i=1:nip %gauss_pts_4
                der=shape_der(points,i,ndim,nod);
                jac=der*coord;
                Det=det(jac);
                deriv=jac\der;
                bee=beemat( deriv,nst );
                eps=bee*eld;
                eps=eps-evpt(:,i,iel);
                sigma=dee*eps;
                stress=sigma+tensor(:,i,iel) ;
                [sigm,dsbar,lode_theta]=invar(stress);
                %  !-----------------------check whether yield is violated-------------------
                f=mocouf(phif,cf,sigm,dsbar,lode_theta);
                if converged==1 && iters==limit
                    devp=stress;
                else
                    if f>=zero
                        [dq1,dq2,dq3]=mocouq(psif,dsbar,lode_theta);
                        [m1,m2,m3]= formm(stress);
                        flow=f*(m1*dq1+m2*dq2+m3*dq3);
                        erate=flow*stress;
                        evp=erate*dt;
                        evpt(:,i,iel)=evpt(:,i,iel)+evp;
                        devp=dee*evp;
                    end
                end
                if f>=0 || (converged==1 || iters==limit)
                    eload=bee'*devp;
                    bload=bload+eload*Det*weights(i);
                end
            end
%             if etype(iel)==2 && f>0
%                 disp('rock is failure');
%             end
            % !-----------------------compute the total bodyloads vector----------------
            for jj=1:length(g)
                if g(jj)~=0
                    bdylds(g(jj))=bdylds(g(jj))+bload(jj);
            % react(g(jj))=react(g(jj))+rload(jj);
                end
            end
%         end
        end %elements_3
        if converged==1 || iters==limit
            break;
        end
    end
    % out put maximum displacement
dmax=max(sqrt(disps(:,1).^2+disps(:,2).^2));
disp(['max displacemnet is ',num2str(dmax)]);
 
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
    fprintf(fid,'%f %f %f\n',(disps(i,1)),(disps(i,2)),0);
end
fprintf(fid,'CELL_DATA %d\n',nels);
fprintf(fid,'SCALARS group double 1\n');
fprintf(fid,'LOOKUP_TABLE table1\n');
for i =1:nels
    fprintf(fid,'%d\n',etype(i));
end
fclose(fid);