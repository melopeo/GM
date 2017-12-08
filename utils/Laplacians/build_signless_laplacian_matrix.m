function Q = build_signless_laplacian_matrix(W, shift)
% L = build_signless_Laplacian_matrix(W, shift)
% This function builds a sparse signless Laplacian Matrix
% input:  W (matrix)    : adjacency matrix
%         shift (scalar): diagonal shift on Laplacian
% output: Q (matrix)    : sparse signless Laplacian matrix
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

if nargin < 2
	shift = 0;
end

n            = size(W,1);              % size of graph
d            = sum(W,2);               % degree vector
dInv         = 1./d;                   % d^{-1}        
dInv(d == 0) = 0;                      % take care of isolated nodes

dInv         = dInv.^(0.5);            % d^{-1/2}
DInv         = spdiags(dInv, 0, n, n); % D^{-1/2}

% Symmetric signless Laplacian
Q = spdiags(d ~= 0, 0, n, n) + DInv*W*DInv; 

% enforce symmetry
Q = (Q  + Q')/2;

if shift > 0
	Q = Q  + shift*speye(n);
end
