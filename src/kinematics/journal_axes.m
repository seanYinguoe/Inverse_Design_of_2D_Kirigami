function [fig, ax] = journal_axes(n, STY, sz, layout)
%JOURNAL_AXES  A white figure with styled, gridless axes in the project style.
%
%   [fig, ax] = journal_axes(n, STY)
%   [fig, ax] = journal_axes(n, STY, [w h])
%   [fig, ax] = journal_axes(n, STY, [w h], [rows cols])
%
%   Shared version of the helper that FORWARD_ANALYSIS defines locally, so that
%   every figure in the forward-analysis family carries identical styling. That
%   script keeps its own local copy - a local function shadows this one - so it
%   is untouched by the existence of this file.
%
%   INPUTS
%     n      : figure number
%     STY    : style struct with fields font, fs, lwAxis (see FORWARD_ANISOTROPY)
%     sz     : [width height] in inches, default [3.5 3.0]
%     layout : [rows cols] of subplots, default [1 1]
%
%   OUTPUT
%     fig : the figure handle
%     ax  : rows-by-cols array of axes handles (scalar for the default layout)

if nargin < 3 || isempty(sz),     sz     = [3.5 3.0]; end
if nargin < 4 || isempty(layout), layout = [1 1];     end

fig = figure(n); clf;
set(fig,'Color','w','Units','inches', ...
    'Position',[1 + 0.35*(n-1), 1, sz(1), sz(2)], ...
    'DefaultAxesFontName',STY.font,'DefaultTextFontName',STY.font, ...
    'DefaultAxesFontSize',STY.fs,'DefaultTextFontSize',STY.fs);

nr = layout(1);   nc = layout(2);
ax = gobjects(nr,nc);
for k = 1:nr*nc
    [j,i] = ind2sub([nc nr], k);        % row-major over the subplot grid
    ax(i,j) = subplot(nr, nc, k, 'Parent', fig);
    hold(ax(i,j),'on');
    set(ax(i,j),'FontName',STY.font,'FontSize',STY.fs,'LineWidth',STY.lwAxis, ...
        'Box','on','TickDir','in','TickLength',[0.018 0.018], ...
        'XColor','k','YColor','k','Layer','top');
    grid(ax(i,j),'off');
end
if isscalar(ax), ax = ax(1); end
end
