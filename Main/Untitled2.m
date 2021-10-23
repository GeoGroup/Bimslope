%DEMO5 assemble triangulations for multi-part geometry defi-
%nitions.
addpath('..\mesh2d');
    fprintf(1, [ ...
' Both the REFINE2 and SMOOTH2 routines also support "multi-  \n', ...
' part" geometry definitions -- generating conforming triang- \n', ...
' ulations that conform to internal and external constraints. \n', ...
        ] ) ;

%---------------------------------------------- create geom.
    nod1 = [
        -1., -1.; +1., -1.
        +1., +1.; -1., +1.
        ] ;
    edg1 = [
         1 ,  2 ;  2 ,  3
         3 ,  4 ;  4 ,  1
        ] ;
    edg1(:,3) = +0;
    
    
    nod2 = [
        +.1, +0.; +.8, +0.
        +.8, +.8; +.1, +.8
        ] ;
    edg2 = [
         1 ,  2 ;  2 ,  3
         3 ,  4 ;  4 ,  1
        ] ;
    edg2(:,3) = +1;
        
        
%     adel = 2.*pi / +64 ;
%     amin = 0.*pi ;
%     amax = 2.*pi - adel;
%     
%     xcir = +.33 * ...
%         cos(amin:adel:amax)';
%     ycir = +.33 * ...
%         sin(amin:adel:amax)';
%     xcir = xcir - .33;
%     ycir = ycir - .25;
%     ncir = [xcir,ycir] ;
%     numc = size(ncir,1);
%         
%     ecir(:,1) = ...
%         [(1:numc-1)'; numc] ;
%     ecir(:,2) = ...
%         [(2:numc-0)'; +1  ] ;
%     ecir(:,3) = +2;
    
    edg2(:,1:2) = ...
    edg2(:,1:2)+size(nod1,1);
    edge = [edg1; edg2];
    node = [nod1; nod2];
        
%     ecir(:,1:2) = ...
%     ecir(:,1:2)+size(node,1);
%     edge = [edge; ecir];
%     node = [node; ncir];
    
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
    hmax = +0.1 ;
 
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
    
%     figure;
%     patch('faces',tlfs(:,1:3),'vertices',vlfs , ...
%         'facevertexcdata' , hlfs, ...
%         'facecolor','interp', ...
%         'edgecolor','none') ;
%     hold on; axis image off;
%     title(['MESH-SIZE: KIND=DELAUNAY, |TRIA|=', ...
%         num2str(size(tlfs,1))]) ;
%         
%     tricost(vert,etri,tria,tnum);
%            
%     drawnow;
%         
%     set(figure(1),'units','normalized', ...
%         'position',[.05,.50,.30,.35]) ;
%     set(figure(2),'units','normalized', ...
%         'position',[.35,.50,.30,.35]) ; 
%     set(figure(3),'units','normalized', ...
%         'position',[.05,.05,.30,.35]) ;