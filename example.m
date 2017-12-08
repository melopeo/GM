function example

addpath(genpath('utils'))

%% Toy example where both graphs are noiseless
n           = 100;
Wpos        = sparse([ones(n), zeros(n); zeros(n), ones(n)]);
Wneg        = sparse([zeros(n), ones(n); ones(n), zeros(n)]);
numClusters = 2;
shift       = 1.e-6;
C           = clustering_signed_networks_with_geometric_mean_of_Laplacians(Wpos,Wneg,numClusters,shift);
figure, stem(C)
1;

%% Toy example where the positive graph is noiseless and negative graph is not informative
n           = 100;
Wpos        = sparse([ones(n), zeros(n); zeros(n), ones(n)]);
Wneg        = sparse(ones(2*n));
numClusters = 2;
shift       = 1.e-6;
C           = clustering_signed_networks_with_geometric_mean_of_Laplacians(Wpos,Wneg,numClusters,shift);
figure, stem(C)
1;

%% Toy example where the positive graph is not informative and negative graph is noiseless
n           = 100;
Wpos        = sparse(ones(2*n));
Wneg        = sparse([zeros(n), ones(n); ones(n), zeros(n)]);
numClusters = 2;
shift       = 1.e-6;
C           = clustering_signed_networks_with_geometric_mean_of_Laplacians(Wpos,Wneg,numClusters,shift);
figure, stem(C)
1;