function L = get_geometric_mean_Laplacian(Wpos,Wneg,shift)
% C = get_geometric_mean_Laplacian(Wpos,Wneg,shift)
% This function builds the geometric mean Laplacian
% Input:      Wpos (matrix)             : adjacency matrix of positive/friendly relationships
%             Wneg (matrix)             : adjacency matrix of negative/enmisty relationships
%             shift (scalar) (optional) : diagonal shift on Laplacians
%                                         (default: shift = 1.e-6)
% Reference:
% Mercado, Pedro, Francesco Tudisco, and Matthias Hein. 
% "Clustering Signed Networks with the Geometric Mean of Laplacians." 
% Advances in Neural Information Processing Systems. 2016.

if nargin < 3
    shift = 1.e-6;
end

Lpos = build_laplacian_matrix(Wpos, shift);          % Get Laplacian
Qneg = build_signless_laplacian_matrix(Wneg, shift); % Get signless Laplacian
L    = sharp(full(Qneg), full(Lpos), 1/2);           % Get Geometric Mean Laplacian
