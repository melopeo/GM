function [legendCell, colorCell, markerCell, LineStyleCell, LineWidthCell, MarkerSizeCell] = get_plot_parameters

%% plotting parameters
legendCell     ={'$L^+_{\textrm{sym}}$', '$Q^-_{\textrm{sym}}$', '$L_{\textrm{SN}}$', '$L_{\textrm{BN}}$', '$L_{\textrm{AM}}$', '$L_{\textrm{GM}}(ours)$', '$L_{\textrm{GM}}(ours)(krylov)$'};

numModels      = length(legendCell);
colorMatrix    = distinguishable_colors(numModels);
colorMatrix    = flipud(colorMatrix);
colorCell      = mat2cell(colorMatrix, ones(numModels,1), size(colorMatrix,2));

colorCell      = {'k', [255,127,0]/255, [152,78,163]/255, [77,175,74]/255, [55,126,184]/255, [228,26,28]/255, 'm'};

markerCell     = {'>', '<', '+', 'p', 'o', 's', '*' };
LineStyleCell  = repmat({'-'},1,numModels);
LineWidthCell  = repmat({1.5},1,numModels);
MarkerSizeCell = repmat({5},1,numModels);

