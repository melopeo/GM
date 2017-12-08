function W = create_sparse_matrix_wikipedia

% Load Data
% file_name = 'wikiElec.ElecBs3_without_header.txt';
% [E,VarName2,VarName3,VarName4,VarName5] = importfile(file_name);
% save('wiki_txt_data.mat');
% load('wiki_txt_data.mat');
% 
% 
% Uloc = find(strcmp(E, 'U'));
% Vloc = find(strcmp(E, 'V'));
% 
% data1 = [];
% data2 = [];
% data3 = [];
% 
% for i = 1:length(Uloc)-1
% 
% % i = 1;
% 
% Vloc_temp    = Vloc( (Vloc > Uloc(i)) & (Vloc < Uloc(i+1)) );
% num_voters   = length(Vloc_temp);
% 
% nominee_id   = VarName2( Uloc(i) );
% voter_id     = VarName3( Vloc_temp );
% vote         = VarName2( Vloc_temp );
% 
% 
% data1 = [voter_id; data1];
% data2 = [nominee_id*ones(num_voters, 1); data2];
% data3 = [vote; data3];
% 
% end
% 1;
% save('wiki_txt_data_formatted.mat');
% load('wiki_txt_data_formatted.mat');
% 
% 
% m = max( [ max(data1) max(data2) ] );
% W         = sparse(data1, data2, data3, m, m);
% 
% save('wiki_matrix.mat', 'W');
load('wiki_matrix.mat', 'W');
1;