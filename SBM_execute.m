function SBM_execute

dirName_Output_Data = 'SBM_execute';
if ~exist(dirName_Output_Data,'dir')
    mkdir(dirName_Output_Data)
end

addpath(genpath('utils'))
addpath(genpath('utilsPlots'))
addpath(genpath('utilsSBM'))


% % % % % % % % % % % % % % experiments % % % % % % % % % % % % % % 

numRuns                        = 50;
numClusters                    = 2;
sizeOfEachCluster              = 100*ones(numClusters,1);
diffArray                      = linspace(0.1, -0.1, 30);
Laplacian_str_cell             = {'Laplacian_positive', 'SignlessLaplacian_negative', 'signed_normalized_cut', 'balance_normalized_cut', 'arithmetic_mean', 'geometric_mean'};

errorPerLaplacianMatrixCellPos = SBM_Wpos_fixed(numRuns, numClusters, sizeOfEachCluster, diffArray, Laplacian_str_cell);
errorPerLaplacianMatrixCellNeg = SBM_Wneg_fixed(numRuns, numClusters, sizeOfEachCluster, diffArray, Laplacian_str_cell);

filename = strcat(dirName_Output_Data, filesep, 'output.mat');
save(filename, 'errorPerLaplacianMatrixCellPos', 'errorPerLaplacianMatrixCellNeg', '-v7.3')
% load(filename)
1;

% % % % % % % % % % % % % % plots % % % % % % % % % % % % % % 
[legendCell, colorCell, markerCell, LineStyleCell, LineWidthCell, MarkerSizeCell] = get_plot_parameters;
fontSize        = 20;
fontSize_legend = 14;

fig_handle = figure;
hold on

% Plot: Experiment fixing W^+
plot_SBM_Wpos_Fixed(errorPerLaplacianMatrixCellPos, diffArray, ...
    fontSize, fontSize_legend, legendCell, colorCell, markerCell, LineStyleCell, LineWidthCell, MarkerSizeCell);

% Plot: Experiment fixing W^-
plot_SBM_Wneg_Fixed(errorPerLaplacianMatrixCellNeg, diffArray, ...
    fontSize, fontSize_legend, legendCell, colorCell, markerCell, LineStyleCell, LineWidthCell, MarkerSizeCell);

filename_plot = strcat(dirName_Output_Data, filesep, 'plot');
save_plots(fig_handle, filename_plot)


function errorPerLaplacianMatrixCell = run_experiment_loop(numRuns, pin_positive_Array, pout_positive_Array, pin_negative_Array, pout_negative_Array, diffArray, Laplacian_str_cell, numClusters, sizeOfEachCluster)

    for k = 1:length(diffArray)

        pin_positive  = pin_positive_Array(k);
        pout_positive = pout_positive_Array(k);
        pin_negative  = pin_negative_Array(k);
        pout_negative = pout_negative_Array(k);

        for j = 1:numRuns
            for i = 1:length(Laplacian_str_cell)

                Laplacian_str            = Laplacian_str_cell{i};
                s                        = RandStream('mcg16807', 'Seed', j); RandStream.setGlobalStream(s);
                [C, GroundTruth]         = run_experiment(sizeOfEachCluster, pin_positive, pout_positive, pin_negative, pout_negative, numClusters, Laplacian_str);

                clustering_error         = classification_error_for_clustering(C, GroundTruth);
                errorPerLaplacian(i)     = clustering_error;

            end
            errorPerLaplacianMatrix(j,:) = errorPerLaplacian;
        end
        errorPerLaplacianMatrixCell{k}   = errorPerLaplacianMatrix;
    end
    1;

1;

function [C, GroundTruth] = run_experiment(sizeOfEachclass, pin_positive, pout_potisive, pin_negative, pout_negative, numClusters, Laplacian_str)

% build ground truth array
numberOfClasses = length(sizeOfEachclass);
GroundTruth = [];
for i = 1:numberOfClasses
    GroundTruth = [GroundTruth; i*ones(sizeOfEachclass(i),1)];
end

% generate random graphs
isConnected = 0;
while ~isConnected
    [Wpos,isConnected] = generate_sbm_graph(GroundTruth', pin_positive, pout_potisive);
end

isConnected = 0;
while ~isConnected
    [Wneg,isConnected] = generate_sbm_graph(GroundTruth', pin_negative, pout_negative);
end

if strcmp(Laplacian_str, 'geometric_mean')
    C = clustering_signed_networks_with_geometric_mean_of_Laplacians(Wpos,Wneg,numClusters);
else
    % get Laplacian
    L = get_signed_Laplacian(Wpos,Wneg,Laplacian_str);

    % get eigenvectors
    [eigvecs, ~] = eigs(L, numClusters, 'sa');

    % get clusters
    C = kmeans(eigvecs, numClusters, 'Replicates', 10, 'emptyaction', 'singleton');    
end
    

1;

function errorPerLaplacianMatrixCell = SBM_Wpos_fixed(numRuns, numClusters, sizeOfEachCluster, diffArray, Laplacian_str_cell)
% Experiment fixing W^+

pin_positive_global         = 0.09;
pout_positive_global        = 0.01;

pin_positive_Array          = pin_positive_global*ones(size(diffArray));
pout_positive_Array         = pout_positive_global*ones(size(diffArray));
pin_negative_Array          = (0.1-diffArray)/2;
pout_negative_Array         = (0.1+diffArray)/2;

errorPerLaplacianMatrixCell = run_experiment_loop(numRuns, pin_positive_Array, pout_positive_Array, pin_negative_Array, pout_negative_Array, diffArray, Laplacian_str_cell, numClusters, sizeOfEachCluster);

function errorPerLaplacianMatrixCell = SBM_Wneg_fixed(numRuns, numClusters, sizeOfEachCluster, diffArray, Laplacian_str_cell)
% Experiment fixing W^-

pin_negative_global         = 0.01;
pout_negative_global        = 0.09;

pin_positive_Array          = (0.1+diffArray)/2;
pout_positive_Array         = (0.1-diffArray)/2;
pin_negative_Array          = pin_negative_global*ones(size(diffArray));
pout_negative_Array         = pout_negative_global*ones(size(diffArray));

errorPerLaplacianMatrixCell = run_experiment_loop(numRuns, pin_positive_Array, pout_positive_Array, pin_negative_Array, pout_negative_Array, diffArray, Laplacian_str_cell, numClusters, sizeOfEachCluster);

function plot_SBM_Wpos_Fixed(errorPerLaplacianMatrixCell, diffArray, ...
    fontSize, fontSize_legend, legendCell, colorCell, markerCell, LineStyleCell, LineWidthCell, MarkerSizeCell)

% Experiment fixing W^+
pin_positive_global         = 0.09;
pout_positive_global        = 0.01;

% get means
error_matrix     = cell2mat(cellfun(@(x) median(x,1), errorPerLaplacianMatrixCell, 'UniformOutput', false)');

legendLocation   = 'northwest';

loc              = (diffArray>-0.06).*(diffArray<0.06) == 1;
diffArray        = diffArray(loc);
error_matrix     = error_matrix(loc,:);

error_matrix     = flipud(error_matrix);
title_str        = ['$', 'p^{+}_{\textrm{in}}=',  num2str(pin_positive_global), ',\,\,', 'p^{+}_{\textrm{out}}=', num2str(pout_positive_global), '$'];
title_str        = {'Positive Informative', title_str};
xAxisTitle_str   = '$p^{-}_{\textrm{in}} - p^{-}_{\textrm{out}}$';
yAxisTitle_str   = 'Median Clustering Error';

subplot(1,2,2),
hold on

plot_curves(diffArray, error_matrix, ...
        colorCell, markerCell, LineWidthCell, LineStyleCell, MarkerSizeCell, ...
        title_str, xAxisTitle_str, yAxisTitle_str, ...
        legendCell, legendLocation, fontSize_legend, fontSize)
    
function plot_SBM_Wneg_Fixed(errorPerLaplacianMatrixCell, diffArray, ...
    fontSize, fontSize_legend, legendCell, colorCell, markerCell, LineStyleCell, LineWidthCell, MarkerSizeCell)

% Experiment fixing W^-
pin_negative_global  = 0.01;
pout_negative_global = 0.09;

% get means
error_matrix = cell2mat(cellfun(@(x) median(x,1), errorPerLaplacianMatrixCell, 'UniformOutput', false)');

legendLocation   = 'northeast';

loc              = (diffArray>-0.06).*(diffArray<0.06) == 1;
diffArray        = diffArray(loc);
error_matrix     = error_matrix(loc,:);

% error_matrixPlot = flipud(error_matrixPlot);
title_str        = ['$', 'p^{-}_{\textrm{in}}=',  num2str(pin_negative_global), ',\,\,', 'p^{-}_{\textrm{out}}=', num2str(pout_negative_global), '$'];
title_str        = {'Negative Informative', title_str};
xAxisTitle_str   = '$p^{+}_{\textrm{in}} - p^{+}_{\textrm{out}}$';
yAxisTitle_str   = 'Median Clustering Error';

subplot(1,2,1),
hold on

plot_curves(diffArray, error_matrix, ...
        colorCell, markerCell, LineWidthCell, LineStyleCell, MarkerSizeCell, ...
        title_str, xAxisTitle_str, yAxisTitle_str, ...
        legendCell, legendLocation, fontSize_legend, fontSize)
    
function plot_curves(diffArray, error_matrix, ...
        colorCell, markerCell, LineWidthCell, LineStyleCell, MarkerSizeCell, ...
        title_str, xAxisTitle_str, yAxisTitle_str, ...
        legendCell, legendLocation, fontSize_legend, fontSize)
    
    for i = 1:size(error_matrix,2)
        plot(diffArray, error_matrix(:,i), ...
        'Color',           colorCell{i}, ...
        'Marker',          markerCell{i}, ...
        'MarkerFaceColor', colorCell{i}, ...
        'MarkerEdgeColor', colorCell{i}, ...
        'LineWidth',       LineWidthCell{i}, ...
        'LineStyle',       LineStyleCell{i}, ...
        'MarkerSize',      MarkerSizeCell{i});
    end
    
    title(title_str, 'interpreter', 'latex', 'FontSize', fontSize)   
    xlabel(xAxisTitle_str, 'interpreter', 'latex', 'FontSize', fontSize,'fontweight','bold')
    ylabel(yAxisTitle_str, 'interpreter', 'latex', 'FontSize', fontSize,'fontweight','bold')
    
    axis square tight
    box on
    daspect
    
    set(gca,'XTick',     [-5 0 5]/100)
    set(gca,'XTickLabel',[-5 0 5]/100)
    set(gca,'YTick',0:0.1:0.5)
    
    ylim([0 0.52])
    
    legend(legendCell,'Location',legendLocation, 'Interpreter','latex','FontSize',fontSize_legend)
    
function save_plots(fig_handle, filename_prefix)

    filename   = strcat(filename_prefix, '.fig');
    savefig(filename)
    file_name_pdf = strcat(filename_prefix, '.pdf');
    file_name_png = strcat(filename_prefix, '.png');
    set(fig_handle,'PaperOrientation','landscape');
    fillPage(fig_handle)
    print(fig_handle, '-dpdf', '-r150', file_name_pdf);
    print(fig_handle, '-dpng', '-r150', file_name_png);