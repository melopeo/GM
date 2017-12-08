function [fig_handle_blue_notZoom, fig_handle_red_notZoom, ...
    fig_handle_blue_Zoom, fig_handle_red_Zoom] = ...
    get_spy_plots(Wpos, Wneg, idxSort, C)

% figuresVisible_str = 'off';
figuresVisible_str = 'on';


% Plot positive graph
fig_handle_blue_notZoom = figure('visible',figuresVisible_str); 
% set(gca,'visible','off')
hold on
spyBlue(Wpos(idxSort, idxSort), 'b.')
set(gca,'XTick',[])
set(gca,'YTick',[])
xlabel('')
box on

% Plot negative graph
fig_handle_red_notZoom = figure('visible',figuresVisible_str); 
% set(gca,'visible','off')
hold on
spyRed(Wneg(idxSort, idxSort), 'r.')
set(gca,'XTick',[])
set(gca,'YTick',[])
xlabel('')
box on

% Plot positive graph with ZOOM
fig_handle_blue_Zoom = figure('visible',figuresVisible_str); 
% set(gca,'visible','off')
hold on
spyBlue(Wpos(idxSort(sum(C==1)+1:end), idxSort(sum(C==1)+1:end)), 'b.')
set(gca,'XTick',[])
set(gca,'YTick',[])
xlabel('')
box on

% Plot negative graph with ZOOM
fig_handle_red_Zoom = figure('visible',figuresVisible_str); 
% set(gca,'visible','off')
hold on
spyRed(Wneg(idxSort(sum(C==1)+1:end), idxSort(sum(C==1)+1:end)), 'r.')
set(gca,'XTick',[])
set(gca,'YTick',[])
xlabel('')
box on
1;