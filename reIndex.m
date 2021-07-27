function edges = reIndex(edges)
vertexName = unique(edges, 'stable');
for i=1:length(edges)
    edges(i, 1) = find(vertexName==edges(i, 1));
    edges(i, 2) = find(vertexName==edges(i, 2));
end
