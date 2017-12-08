function L = get_balance_normalized_Laplacian(Wpos,Wneg,shift)
% C = get_balance_normalized_Laplacian(Wpos,Wneg,shift)
% This function builds the balance normalized Laplacian.
% Input:      Wpos (matrix)             : adjacency matrix of positive/friendly relationships
%             Wneg (matrix)             : adjacency matrix of negative/enmisty relationships
%             shift (scalar) (optional) : diagonal shift on Laplacians
%                                         (default: shift = 0)
% Reference:
% Chiang, Kai-Yang, Joyce Jiyoung Whang, and Inderjit S. Dhillon. 
% "Scalable clustering of signed networks using balance normalized cut." 
% Proceedings of the 21st ACM international conference on Information and 
% knowledge management. ACM, 2012.

% Observation: in the paper they assume that edges are either positve or
% negative, i.e. the same edge can not appear in both Wpos and Wneg. Thats
% why we build an adjacency matrix as W = Wpos - Wneg and take the absolute
% degree correspondingly.


if nargin < 3
    shift = 0;
end

W            = Wpos - Wneg;
d            = sum( W.*(W>0),2 ) + sum( abs(W.*(W<0)),2 );
dPos         = sum( W.*(W>0),2 );

n            = size(W,1);                   % size of graph
dInv         = 1./d;                        % d^{-1}                   
dInv(d == 0) = 0;                           % take care of isolated nodes

dInv         = dInv.^(0.5);                 % d^{-1/2}
DInv         = spdiags(dInv, 0, n, n);      % D^{-1/2}
Dpos         = spdiags(dPos, 0, n, n);      % Dpos


L            = DInv*(Dpos - W)*DInv;        % Symmetric Laplacian
L            = (L  + L')/2;                 % enforce symmetry

% add shift
if shift > 0
	L = L  + shift*speye(n);
end
