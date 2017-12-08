function kernelMatrix = get_connected_bipartite_components(W)

numNodes = size(W,1);
d        = sum(W,1);

[comp,connected,sizes] = connectedComponents(W);
numComponents          = length(sizes);

% check if components are bipartite
compReduced_relabeled_unique        = unique(comp);
compReduced_relabeled_unique_length = length(compReduced_relabeled_unique);

isbipartiteArray = zeros(compReduced_relabeled_unique_length,1);
for j = 1:compReduced_relabeled_unique_length
    
    loc           = comp == compReduced_relabeled_unique(j);
    Wtemp         = W(loc, loc);
    adjacencyList = adj2adjL(Wtemp);
    
    [isbipartiteArray(j), A{j}, B{j}]=isbipartite(adjacencyList);
end

% Build Indicator vectors for signless laplacian
numBipartiteComponents = sum(isbipartiteArray);
kernelMatrix           = zeros(numNodes,numBipartiteComponents);
bipartiteIdx           = find(isbipartiteArray);

for j = 1:numBipartiteComponents
    idx = find(comp == compReduced_relabeled_unique(bipartiteIdx(j)));
    v = zeros(numNodes, 1);
    v(idx(A{bipartiteIdx(j)})) = 1;
    v(idx(B{bipartiteIdx(j)})) = -1;
    kernelMatrix(:,j) = v;
end
1;
