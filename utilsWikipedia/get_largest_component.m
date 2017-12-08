function [W_out, loc, connected] = get_largest_component(W)
% W_out = get_largest_component(W)
% This function gets the largest connected component of a given adjacency
% matrix
% INPUT : W     : adjacency matrix
% OUTPUT: W_out : adjacency matrix of largest connected component
[comp,connected,sizes] = connectedComponents(abs(W));
idx                    = find(sizes == max(sizes), 1, 'first');
loc                    = comp == idx;
W_out                  = W(loc, loc);

if ~connected
    warning('get_largest_component: Input graph is not connected');
end