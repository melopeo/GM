function L = get_signed_normalized_Laplacian(Wpos,Wneg,shift)
% C = get_signed_normalized_cut(Wpos,Wneg,shift)
% This function builds the signed Laplacian.
% Input:      Wpos (matrix)             : adjacency matrix of positive/friendly relationships
%             Wneg (matrix)             : adjacency matrix of negative/enmisty relationships
%             shift (scalar) (optional) : diagonal shift on Laplacians
%                                         (default: shift = 0)
% Reference:
% Kunegis, Jérôme, et al. 
% "Spectral analysis of signed graphs for clustering, prediction and visualization." 
% Proceedings of the 2010 SIAM International Conference on Data Mining. 
% Society for Industrial and Applied Mathematics, 2010.

% Observation: in the paper they assume that edges are either positve or
% negative, i.e. the same edge can not appear in both Wpos and Wneg. Thats
% why we build an adjacency matrix as W = Wpos - Wneg and take the absolute
% degree correspondingly.


if nargin < 3
    shift = 0;
end

W            = Wpos - Wneg;
d            = sum( W.*(W>0),2 ) + sum( abs(W.*(W<0)),2 );

n            = size(W,1);                   % size of graph
dInv         = 1./d;                        % d^{-1}                   
dInv(d == 0) = 0;                           % take care of isolated nodes

dInv         = dInv.^(0.5);                 % d^{-1/2}
DInv         = spdiags(dInv, 0, n, n);      % D^{-1/2}

L = spdiags(d ~= 0, 0, n, n) - DInv*W*DInv; % Symmetric Laplacian
L = (L  + L')/2;                            % enforce symmetry

% add shift
if shift > 0
	L = L  + shift*speye(n);
end
