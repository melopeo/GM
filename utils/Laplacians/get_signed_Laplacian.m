function L = get_signed_Laplacian(Wpos,Wneg,Laplacian_str,shift)
% C = get_arithmetic_mean_Laplacian(Wpos,Wneg,shift)
% This function builds the arighmetic mean Laplacian for signed networks
% Input:      Wpos (matrix)             : adjacency matrix of positive/friendly relationships
%             Wneg (matrix)             : adjacency matrix of negative/enmisty relationships
%             Laplacian_str (string)    : kind of signed Laplacian to build
%             shift (scalar) (optional) : diagonal shift on Laplacians
%                                         (default: shift = 0)

if nargin < 4
    shift = 0;
end

if nargin < 4 
    if strcmp(Laplacian_str, 'geometric_mean')
        shift = 1.e-6;
    else
        shift = 0;
    end 
end

switch Laplacian_str
    
    case 'arithmetic_mean'
        if nargin < 4, shift = 0;     end
        L = get_arithmetic_mean_Laplacian(Wpos,Wneg,shift);
        
%     case 'geometric_mean'
%         if nargin < 4, shift = 1.e-6; end
%         L = get_geometric_mean_Laplacian(Wpos,Wneg,shift);
        
    case 'signed_normalized_cut'
        if nargin < 4, shift = 0;     end
        L = get_signed_normalized_Laplacian(Wpos,Wneg,shift);
        
    case 'balance_normalized_cut'
        if nargin < 4, shift = 0;     end
        L = get_balance_normalized_Laplacian(Wpos,Wneg,shift);
        
    case 'Laplacian_positive'
        if nargin < 4, shift = 0;     end
        L = build_laplacian_matrix(Wpos, shift);
        
    case 'SignlessLaplacian_negative'
        if nargin < 4, shift = 0;     end
        L = build_signless_laplacian_matrix(Wneg, shift);
        
    otherwise
        warning('Unexpected Laplacian_str value.')

end

