function collision_particle=poly_collistion(poly,polygon,particlenum)
% Judge the intersection of a polygon and a set of polygons
collision_particle=0;
poly=polyshape(poly(:,1),poly(:,2));
if particlenum==0
    collision_particle=0;
else
    for i=1:particlenum
        polyout = separation_axis_theorem_2d(get_convhull_2d(poly.Vertices), get_convhull_2d(polygon{i}.Vertices)  );
        if polyout 
            collision_particle=1;
            break;
        end
    end
end




end

