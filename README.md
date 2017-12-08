## Clustering Signed Networks with the Geometric Mean of Laplacians

MATLAB implementation of the paper:

[P. Mercado, F. Tudisco, and M. Hein, Clustering signed networks with the geometric mean of Laplacians. In NIPS 2016.](http://papers.nips.cc/paper/6164-clustering-signed-networks-with-the-geometric-mean-of-laplacians.pdf)

## Content:
- `example.m` : contains three easy examples showing how to use the code

- `SBM_execute.m` : runs spectral clustering on signed networks under the stochastic block model (Fig.2 of our [paper](http://papers.nips.cc/paper/6164-clustering-signed-networks-with-the-geometric-mean-of-laplacians.pdf))

- `wikipedia_example_analysis.m` : generates plots from Wikipedia experiment (Fig.4 of our [paper](http://papers.nips.cc/paper/6164-clustering-signed-networks-with-the-geometric-mean-of-laplacians.pdf)). 
   - The corresponding eigenvectors are precomputed with `wikipedia_example.m` and are located in the folder `Wikipedia/wikipedia_example`
   
## Usage:
Let `Wpos` and `Wneg` be the positive and negative adjacency matrices, respectively, and `numClusters` the desired number of clusters. Clusters through the geometric mean Laplacian are computed via

```
C = clustering_signed_networks_with_geometric_mean_of_Laplacians(Wpos,Wneg,numClusters);
```
## Quick Overview
![](https://github.com/melopeo/GM/blob/master/PaperAndPoster/ClusteringSignedNetworksWithTheGeometricMeanOfLaplaciansPoster.jpg)

# Citation:
```
@incollection{NIPS2016_6164,
title = {Clustering Signed Networks with the Geometric Mean of Laplacians},
author = {Mercado, Pedro and Tudisco, Francesco and Hein, Matthias},
booktitle = {Advances in Neural Information Processing Systems 29},
editor = {D. D. Lee and M. Sugiyama and U. V. Luxburg and I. Guyon and R. Garnett},
pages = {4421--4429},
year = {2016},
publisher = {Curran Associates, Inc.},
url = {http://papers.nips.cc/paper/6164-clustering-signed-networks-with-the-geometric-mean-of-laplacians.pdf}
}
```
