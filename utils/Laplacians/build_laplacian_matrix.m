function L = build_laplacian_matrix(W, shift)
% L = build_Laplacian_matrix(W, shift)
% This function builds a sparse Laplacian Matrix
% input:  W (matrix)    : adjacency matrix
%         shift (scalar): diagonal shift on Laplacian
% output: L (matrix)    : sparse Laplacian matrix
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

% Symmetric Laplacian
L = spdiags(d ~= 0, 0, n, n) - DInv*W*DInv; 

% enforce symmetry
L = (L  + L')/2;

% add shift
if shift > 0
	L = L  + shift*speye(n);
end
