% FIG_ROTATING_UNITS  Illustration: why squares and triangles tilt and hexagons do not.
%
%   STANDALONE FIGURE SCRIPT - nothing else in the project depends on it, so
%   it can be deleted freely. Run it from this folder.
%
%   THE ARGUMENT
%     In a rotating-units mechanism the rigid panels are pinned to each other
%     at their corners. A pin survives the motion only if the two panels it
%     joins turn by the SAME angle in OPPOSITE senses, so neighbouring panels
%     must carry opposite signs: + - + - ... The panel adjacency graph has to
%     be two-colourable.
%
%     Walk once around a vertex of a regular tiling and that walk is a closed
%     cycle of length z, the number of panels meeting there. A cycle can be
%     two-coloured only if its length is EVEN, hence
%
%       triangles  z = 6  even  ->  + - + - + -   closes      mechanism
%       squares    z = 4  even  ->  + - + -       closes      mechanism
%       hexagons   z = 3  odd   ->  + - ?         contradicts locked
%
%     The third hexagon touches one panel that forced it to be + and another
%     that forced it to be -, so no consistent assignment exists and the
%     honeycomb admits no rotating-units mode. Same reason a 3-cycle cannot
%     be 2-coloured.
%
%   Row 1 draws that vertex count; row 2 draws the consequence, using
%   FORWARD_DEPLOYMENT for the two tilings that do work.
%
%   See also FORWARD_MAIN, FORWARD_DEPLOYMENT, PLOT_FORWARD, FORWARD_COLORS.

clear;
clc;
addpath(fileparts(fileparts(fileparts(mfilename('fullpath'))))); setup_project;

C     = forward_colors();
cPlus = C(1,:);              % (+) counter-clockwise, blue
cMin  = C(2,:);              % (-) clockwise, orange
cBad  = [0.85 0.10 0.10];    % the contradiction
face  = [1 .98 .54];         % same panel fill as PLOT_FORWARD
dead  = [0.86 0.86 0.86];    % a panel that cannot move
a     = 1;

fig = figure(3); clf;
set(fig,'Color','w','Position',[80 80 1240 760]);

%% ----------------------------------------------------- row 1: the vertex
spec = { 3, 6, '6 triangles at a vertex', 'z = 6 even  \rightarrow  + - + - + -  closes'
         4, 4, '4 squares at a vertex',   'z = 4 even  \rightarrow  + - + -  closes'
         6, 3, '3 hexagons at a vertex',  'z = 3 odd  \rightarrow  + - ?  contradicts' };

for s = 1:3
    p = spec{s,1};   z = spec{s,2};
    ax = subplot(2,3,s); hold(ax,'on');
    rin = 0.60*a/(2*tan(pi/p));            % arrow radius, from the inradius

    for k = 0:z-1
        V = poly_at_vertex(a, p, k*2*pi/z);
        c = mean(V,1);
        conflict = (z == 3) && (k == 2);   % the panel the parity argument breaks on
        if conflict
            patch('XData',V(:,1),'YData',V(:,2),'FaceColor',dead, ...
                  'EdgeColor',[.3 .3 .3],'LineWidth',1.2,'Parent',ax);
            % both signs are demanded of it - draw them faintly, then cross out
            draw_spin(ax, c, rin, +1, 1-0.30*(1-cPlus), 1.6);
            draw_spin(ax, c, rin, -1, 1-0.30*(1-cMin ), 1.6);
            d = 0.55*rin;
            plot(ax, c(1)+[-d d], c(2)+[-d d], '-', 'Color',cBad,'LineWidth',4);
            plot(ax, c(1)+[-d d], c(2)+[ d -d], '-', 'Color',cBad,'LineWidth',4);
        else
            patch('XData',V(:,1),'YData',V(:,2),'FaceColor',face, ...
                  'EdgeColor',[.3 .3 .3],'LineWidth',1.2,'Parent',ax);
            if mod(k,2) == 0
                draw_spin(ax, c, rin, +1, cPlus, 2.6);
            else
                draw_spin(ax, c, rin, -1, cMin, 2.6);
            end
        end
    end

    plot(ax, 0, 0, 'o', 'MarkerFaceColor',[.85 .1 .1], ...
        'MarkerEdgeColor','w','LineWidth',1,'MarkerSize',11);
    axis(ax,'equal'); axis(ax,'off');
    title(ax, {spec{s,3}, spec{s,4}}, 'FontSize',11);
end

%% ------------------------------------------------ row 2: the consequence
ax = subplot(2,3,4);
St = forward_deployment('triangle',3,4,a,80*pi/180);
plot_forward(St, 'Parent',ax);
title(ax, {'deploys', sprintf('\\theta = 80\\circ,  area \\times %.1f', St.lambda^2)}, 'FontSize',11);

ax = subplot(2,3,5);
Ss = forward_deployment('square',2,2,a,60*pi/180);
plot_forward(Ss, 'Parent',ax);
title(ax, {'deploys', sprintf('\\theta = 60\\circ,  area \\times %.1f', Ss.lambda^2)}, 'FontSize',11);

ax = subplot(2,3,6); hold(ax,'on');
H = honeycomb_patch(a, 3, 3);
for k = 1:numel(H)
    patch('XData',H{k}(:,1),'YData',H{k}(:,2),'FaceColor',dead, ...
          'EdgeColor',[.35 .35 .35],'LineWidth',1.0,'Parent',ax);
end
axis(ax,'equal'); axis(ax,'off');
title(ax, {'locked', 'no alternating \pm mode'}, 'FontSize',11,'Color',cBad);

sgtitle({'A rotating-units mechanism needs neighbouring panels to spin in opposite senses,', ...
         'so it exists only where an EVEN number of panels meets at every vertex'}, ...
        'FontSize',13,'FontWeight','bold');

% colour key for the two rotation senses
annotation(fig,'textbox',[0.005 0.50 0.99 0.045],'String', ...
    '\color[rgb]{0,0.45,0.74}\bf+\phi  counter-clockwise          \color[rgb]{0.93,0.53,0}\bf-\phi  clockwise          \color[rgb]{0.85,0.1,0.1}\bf\bullet  shared vertex', ...
    'HorizontalAlignment','center','EdgeColor','none','FontSize',11,'Interpreter','tex');

%% ---------------------------------------------------------------- helpers
function V = poly_at_vertex(a, p, phi0)
%POLY_AT_VERTEX  Regular p-gon of side a with one vertex at the origin.
%   Its two edges at the origin leave along phi0 and phi0 + 2*pi/p*(p-2)/2,
%   i.e. it fills the sector of interior angle (p-2)*pi/p starting at phi0.
d = phi0 + (0:p-1)*(2*pi/p);
V = [0 0; cumsum(a*[cos(d(:)) sin(d(:))])];
V = V(1:p,:);
end

function draw_spin(ax, c, r, sense, col, lw)
%DRAW_SPIN  Curved arrow round c, counter-clockwise for sense +1.
sweep = 280;
t  = linspace(90 - sense*sweep/2, 90 + sense*sweep/2, 80)*pi/180;
x  = c(1) + r*cos(t);   y = c(2) + r*sin(t);
plot(ax, x, y, '-', 'Color',col, 'LineWidth',lw);
d  = [x(end)-x(end-1), y(end)-y(end-1)];   d = d/norm(d);
nn = [-d(2) d(1)];
h  = 0.55*r;   w = 0.28*r;
tip = [x(end) y(end)] + d*0.5*h;
patch('XData',[tip(1) tip(1)-h*d(1)+w*nn(1) tip(1)-h*d(1)-w*nn(1)], ...
      'YData',[tip(2) tip(2)-h*d(2)+w*nn(2) tip(2)-h*d(2)-w*nn(2)], ...
      'FaceColor',col,'EdgeColor','none','Parent',ax);
end

function H = honeycomb_patch(a, m, n)
%HONEYCOMB_PATCH  m-by-n patch of regular hexagons of side a (pointy top).
ang = (90 + (0:5)*60)*pi/180;
hex = a*[cos(ang(:)) sin(ang(:))];
H = cell(m*n,1);   k = 0;
for q = 0:m-1
    for p = (0:n-1) - round(q/2)
        c = p*[a*sqrt(3) 0] + q*[a*sqrt(3)/2 1.5*a];
        k = k + 1;
        H{k} = hex + c;
    end
end
end
