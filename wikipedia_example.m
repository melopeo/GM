function wikipedia_example

dirName_Output_Data = strcat('Wikipedia', filesep, 'wikipedia_example');
if ~exist(dirName_Output_Data,'dir')
    mkdir(dirName_Output_Data)
end

addpath(genpath('utilsWikipedia'))
addpath(genpath('utils'))

% load wikipedia matrix
W = create_sparse_matrix_wikipedia;

% We want a symmetric adjacency matrix
W = sign(W + W');

% We take the largest connected component of the positive graph
[~,loc] = get_largest_component( W.*(W>0) );
W       = W(loc, loc);
Wpos    =  W.*(W>0); % adjacency matrix of positive graph
Wneg    = -W.*(W<0); % adjacency matrix of negative graph
 
% Avoid self loops
Wpos = triu(Wpos,1) + triu(Wpos,1)';
Wneg = triu(Wneg,1) + triu(Wneg,1)';

% Get kernels
kernelMatrix_Laplacian         = get_connected_components(Wpos); % Connected components in positive graph
kernelMatrix_singlessLaplacian = get_connected_bipartite_components(Wneg);% Bipartite connected components in negative graph
kernelSize                     = size(kernelMatrix_Laplacian,2) + ... % Global number of components and biparte components
    size(kernelMatrix_singlessLaplacian,2);

% We look for 30 clusters
numEigenvectors             = 30;
numClusters                 = 30;
numEigenvectorsToCompute    = kernelSize+numEigenvectors;

%% Geometric Mean:  L_sym # Q_sym

% Parameters
krylovOpts             = struct;
eigsOpts               = struct;
[krylovOpts, eigsOpts] = process_inputs(krylovOpts, eigsOpts);
parallelExecution      = true;

% Laplacians
shift        = 1.e-6;
L            = build_laplacian_matrix(Wpos, shift);
Q            = build_signless_laplacian_matrix(Wneg, shift);

% Input kernel
kernelMatrix = sum(Wpos, 2).^0.5;

% Get eigenvectors
[eigvecs, eigvals] = power_method_for_geometric_mean_kpik(L, Q, numEigenvectorsToCompute, kernelMatrix, krylovOpts, eigsOpts, parallelExecution);

% sort eigenvalues and eigenvectors
[eigvalsSorted, idxSorted] = sort(eigvals, 'ascend');
eigvecsSorted              = eigvecs(:, idxSorted);

% Take 30 eigenvectors on top of components
eigenvectors                = eigvecsSorted(:,kernelSize+1:kernelSize+numEigenvectors);
1;

filename = strcat(dirName_Output_Data, filesep, 'output2.mat');
save(filename, 'eigenvectors', 'Wpos', 'Wneg', '-v7.3')
 