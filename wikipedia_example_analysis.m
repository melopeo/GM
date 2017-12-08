function wikipedia_example_analysis

dirName_Input_Data  = strcat('Wikipedia', filesep, 'wikipedia_example');
dirName_Output_Data = strcat('Wikipedia', filesep, 'wikipedia_example_analysis');
if ~exist(dirName_Output_Data,'dir')
    mkdir(dirName_Output_Data)
end

addpath(genpath('utilsWikipedia'))
addpath(genpath('utilsPlots'))

% We look for 30 clusters
numClusters = 30;

% load data
filename   = strcat(dirName_Input_Data, filesep, 'output2.mat');
data       = load(filename);

eigenvectors = data.eigenvectors;
Wpos         = data.Wpos;
Wneg         = data.Wneg;

% Apply k-means
randomSeed = 0;
s = RandStream('mcg16807','Seed',randomSeed); RandStream.setGlobalStream(s);
C = kmeans(eigenvectors, numClusters, 'Replicates', 10, 'emptyaction', 'singleton');

% Relabel clusters according to cardinality
C            = relabel_cluster_by_cardinality(C);
[C, idxSort] = sort(C, 'ascend');

[fig_handle_blue_notZoom, fig_handle_red_notZoom, ...
fig_handle_blue_Zoom, fig_handle_red_Zoom] = ...
get_spy_plots(Wpos, Wneg, idxSort, C);
1;
filename_prefix = strcat(dirName_Output_Data, filesep, 'blue_notZoom');
save_plots(fig_handle_blue_notZoom, filename_prefix)

filename_prefix = strcat(dirName_Output_Data, filesep, 'red_notZoom');
save_plots(fig_handle_red_notZoom,  filename_prefix)

filename_prefix = strcat(dirName_Output_Data, filesep, 'blue_Zoom');
save_plots(fig_handle_blue_Zoom,    filename_prefix)

filename_prefix = strcat(dirName_Output_Data, filesep, 'red_Zoom');
save_plots(fig_handle_red_Zoom,     filename_prefix)

