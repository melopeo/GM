function kernelMatrix = get_connected_components_V2(W)

numNodes = size(W,1);

[comp,connected,sizes] = connectedComponents(W);
numComponents          = length(sizes);

% Build Indicator vectors for laplacian
kernelMatrix = zeros(numNodes, numComponents);
for j = 1:numComponents
    loc = comp == j;
    kernelMatrix(:,j) = loc;
end

