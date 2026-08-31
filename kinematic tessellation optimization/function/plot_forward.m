function plot_forward(S, varargin)
%PLOT_FORWARD  Draw one configuration of a forward-deployed tessellation.
%
%   plot_forward(S)
%   plot_forward(S, 'Name', value, ...)
%
%   Draws the rigid panels as filled patches, every shared (cut) edge in the
%   colour it carries for the whole motion, and every hinge as a red dot. The
%   two halves of a cut are drawn separately, so at theta = 0 they sit on top
%   of each other and as theta grows they peel apart into the walls of a void
%   - which is the whole point of colouring them.
%
%   INPUT
%     S : configuration struct from FORWARD_DEPLOYMENT
%
%   OPTIONS
%     'PanelColor'  fill colour of the rigid panels   [1 .98 .54]
%     'EdgeWidth'   line width of the coloured cuts   scaled to the patch
%     'HingeSize'   marker size of the hinge dots     scaled to the patch
%     'ShowHinges'  draw the hinge dots               true
%     'ShowCuts'    draw the coloured cuts            true
%     'Parent'      axes to draw into                 gca
%
%   See also FORWARD_DEPLOYMENT, FORWARD_COLORS, FORWARD_MAIN.

p = inputParser;
p.addParameter('PanelColor', [1 .98 .54]);
p.addParameter('EdgeWidth',  []);
p.addParameter('HingeSize',  []);
p.addParameter('ShowHinges', true);
p.addParameter('ShowCuts',   true);
p.addParameter('Parent',     []);
p.parse(varargin{:});
o  = p.Results;
ax = o.Parent;
if isempty(ax), ax = gca; end
% weight the cuts and the hinge dots by how much of the view one panel spans,
% so a 2-by-2 patch and a 12-by-12 one both stay readable
pf = S.scale/S.view_span;
if isempty(o.EdgeWidth), o.EdgeWidth = min(3.0, max(0.7, 20*pf)); end
if isempty(o.HingeSize), o.HingeSize = min(9.0, max(2.5, 50*pf)); end

cmap = forward_colors();
hold(ax,'on');

% rigid panels
for k = 1:numel(S.panels)
    v = S.panels{k};
    patch('XData',v(:,1),'YData',v(:,2), ...
          'FaceColor',o.PanelColor,'EdgeColor',[.70 .70 .70], ...
          'LineWidth',0.4,'Parent',ax);
end

% shared edges, one line per panel-side of the cut
if o.ShowCuts
    for c = 1:S.ncolors
        k = find(S.edges.cid == c);
        if isempty(k), continue; end
        % NaN-separated so the whole colour goes in as a single line object
        X = [S.edges.seg(k,1) S.edges.seg(k,3) nan(numel(k),1)]';
        Y = [S.edges.seg(k,2) S.edges.seg(k,4) nan(numel(k),1)]';
        plot(ax, X(:), Y(:), '-', 'Color', cmap(c,:), ...
             'LineWidth', o.EdgeWidth, 'DisplayName', sprintf('cut %d', c));
    end
end

% hinges
if o.ShowHinges && ~isempty(S.hinges)
    plot(ax, S.hinges(:,1), S.hinges(:,2), 'o', ...
        'MarkerFaceColor',[.85 .10 .10], 'MarkerEdgeColor','none', ...
        'MarkerSize', o.HingeSize, 'DisplayName','hinge');
end

axis(ax,'equal');
axis(ax,'off');
end
