function C = clustering_signed_networks_with_geometric_mean_of_Laplacians(Wpos,Wneg,numClusters,shift,kernelMatrix,krylovOpts,eigsOpts, parallelExecution)
% C = clustering_signed_networks_with_geometric_mean_of_Laplacians(Wpos,Wneg,numClusters,shift)
% This function performs clustering signed networwks with the geometric
% mean of Laplacians
% Input:      Wpos (matrix)              : adjacency matrix of positive/friendly relationships
%             Wneg (matrix)              : adjacency matrix of negative/enmisty relationships
%             numClusters (scalar)       : number of clusters to identify
%             shift (scalar) (optional)  : diagonal shift on Laplacians
%                                         (default: shift = 1.e-6)
%             kernelMatrix (matrix)      : matrix of precomputed eignevectors
%                                         (default: kernelMatrix=[])
%             krylovOpts (structure)     : structure with parameters of
%                                         Krylov method 
%                                          (default:
%                                          -krylovOpts.solverMaxIter = 200;
%                                          -krylovOpts.numIter = 5;
%                                          -krylovOpts.TOL = 1.e-5;)
%             eigsOpts (structure)       : structure with parameters of
%                                         power method
%                                           (default:
%                                          -eigsOpts.numIter = 1000;
%                                          -krylovOpts.numIter = 5;
%                                          -eigsOpts.TOL = 1.e-6;)
%             parallelExecution (boolean): if true: runs power method in parallel
%                                          
% 
% addpath('powerMethod/')
% addpath('utils/')

assert( issparse(Wpos) && issparse(Wneg), 'Wpos and Wneg must be sparse matrices')

if nargin < 4
    shift = 1.e-6;
end

if nargin < 5
    kernelMatrix = [];
end

if nargin < 6
    krylovOpts = struct;
end

if nargin < 7
    eigsOpts = struct;
end

if nargin < 8
    parallelExecution = false;
end

[krylovOpts, eigsOpts] = process_inputs(krylovOpts, eigsOpts);

L                      = build_laplacian_matrix(Wpos, shift);
Q                      = build_signless_laplacian_matrix(Wneg, shift);

[eigvecs, eigvals]     = power_method_for_geometric_mean_kpik(L, Q, numClusters, kernelMatrix, krylovOpts, eigsOpts, parallelExecution);

C                      = kmeans(eigvecs, numClusters, 'Replicates', 10, 'emptyaction', 'singleton');
