function L = get_arithmetic_mean_Laplacian(Wpos,Wneg,shift)
% C = get_arithmetic_mean_Laplacian(Wpos,Wneg,shift)
% This function builds the arighmetic mean Laplacian for signed networks
% Input:      Wpos (matrix)             : adjacency matrix of positive/friendly relationships
%             Wneg (matrix)             : adjacency matrix of negative/enmisty relationships
%             shift (scalar) (optional) : diagonal shift on Laplacians
%                                         (default: shift = 0)

if nargin < 3
    shift = 0;
end

Lpos    = build_laplacian_matrix(Wpos, shift);          % Get Laplacian
Qneg    = build_signless_laplacian_matrix(Wneg, shift); % Get signless Laplacian
L       = 0.5*(Lpos+Qneg);
