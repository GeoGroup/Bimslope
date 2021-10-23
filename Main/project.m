function min_max = project(vertices, axis)
    projections = zeros(1, length(vertices));
    for i = 1:length(vertices)
        projections(i) = dot(vertices(i, :), axis);
    end
    min_max = [min(projections), max(projections)];
end