function test = separation_axis_theorem_2d(vertices_a, vertices_b)
    dirs_a = vertices_to_edges(vertices_a);
    dirs_b = vertices_to_edges(vertices_b);
    dirs = [dirs_a; dirs_b];
    for i = 1:length(dirs)
        dir = dirs(i, :);
        axis = [dir(2), -dir(1)];
        axis = axis/norm(axis);
        a = project(vertices_a, axis);
        b = project(vertices_b, axis);
        overlapping = overlap(a, b);
        if ~overlapping
            test = false;
            return;
        end
    end
    test = true;
end


