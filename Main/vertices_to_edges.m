function dirs = vertices_to_edges(vertices)
    directions = zeros(size(vertices));
    for i = 1:length(vertices) - 1
        p0 = vertices(i, :);
        p1 = vertices(i+1, :);
        directions(i,:) = (p1 - p0);
    end
    p0 = vertices(length(vertices), :);
    p1 = vertices(1, :);
    directions(length(vertices), :) = (p1 - p0);
    dirs = directions;
end