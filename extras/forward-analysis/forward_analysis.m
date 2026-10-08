% FORWARD_ANALYSIS  Porosity and rotational energy against macroscopic strain.
%
%   SCRIPT. Draws two separate journal-scale figures - figure 1 porosity,
%   figure 2 normalised rotational energy - both against macroscopic strain.
%   Compares the two rotating-units patterns built by
%   FORWARD_DEPLOYMENT - rotating squares and rotating triangles - along
%   their whole deployment path, as continuum quantities rather than as
%   pictures.
%
%   MODEL
%     Both patterns have one degree of freedom, the opening angle theta,
%     which is the relative rotation across every hinge. For a tiling of
%     regular polygons of interior angle phi,
%
%       lambda = sin((phi + theta)/2) / sin(phi/2)     linear expansion
%       eps    = lambda - 1                            macroscopic strain
%       p      = 1 - lambda^(-2)                       porosity (void fraction)
%       u      = c * K * theta^2 / l^2                 normalised rotational
%                                                      energy per unit area
%
%     The expansion is isotropic, so the areal expansion is lambda^2 and the
%     Poisson's ratio is -1. The energy is that of a linear torsion spring of
%     stiffness K at every hinge, with c collecting the hinge count per unit
%     cell; theta is in RADIANS, and it is normalised by setting K/l^2 = 1.
%
%       squares    phi =  90 deg,  theta in [0,  90] deg,  c = 1
%       triangles  phi =  60 deg,  theta in [0, 120] deg,  c = sqrt(3)
%
%     The lambda above is the same one SQUARE_DEPLOYMENT and
%     TRIANGLE_DEPLOYMENT use (there it reads cos(theta/2)+sin(theta/2) and
%     2*sin(theta/2+pi/6) respectively); the block at the end of this script
%     checks the two agree.
%
%   See also FORWARD_MAIN, FORWARD_DEPLOYMENT, FORWARD_COLORS.

clear;
clc;
addpath(fileparts(fileparts(fileparts(mfilename('fullpath'))))); setup_project;

K_over_l2 = 1;      % hinge stiffness / (length scale)^2

%% --------------------------------------------------------- the two patterns
pat(1) = struct('name','rotating squares',  'phi',90, 'theta_max',90,  'c',1);
pat(2) = struct('name','rotating triangles','phi',60, 'theta_max',120, 'c',sqrt(3));

C   = forward_colors();
col = C([1 2],:);   % blue for squares, orange for triangles

for k = 1:numel(pat)
    th  = linspace(0, pat(k).theta_max, 601)*pi/180;   % radians
    phi = pat(k).phi*pi/180;

    pat(k).theta  = th;
    pat(k).lambda = sin((phi + th)/2) / sin(phi/2);
    pat(k).eps    = pat(k).lambda - 1;
    pat(k).p      = 1 - pat(k).lambda.^(-2);
    pat(k).u      = pat(k).c * K_over_l2 * th.^2;
end

%% ------------------------------------------------------------ figure style
% One style table for both figures: a single font, a single line weight, no
% grid. Sizes are in points and the figures in inches, so each comes out at
% journal single-column scale without rescaling on export.
STY = struct('font','Helvetica', 'fs',9, 'fsLab',10, 'lw',1.6, ...
             'lwAxis',0.9, 'ms',5);

SAVE_FIG   = false;    % write a vector PDF per figure, next to this script
LABEL_ENDS = false;    % print each endpoint's coordinates next to its marker.
                       % Off by default: the numbers are already listed in the
                       % console table below, and on a journal figure they
                       % belong in the caption rather than over the curve.

% Every curve carries the SAME line weight, so the two porosity curves - which
% are not merely close but identical, p depending on lambda alone - would hide
% each other. They are separated by line style instead: squares solid,
% triangles dashed and drawn second, so the dashes ride visibly along the
% solid blue and the coincidence is the thing you see.
ls    = {'-','--'};
short = {'Squares','Triangles'};
xmax  = 1.05 + 0.22*LABEL_ENDS;

%% ------------------------------------------------------- figure 1: porosity
[fig1, ax1] = journal_axes(1, STY);

h1 = gobjects(numel(pat),1);
for k = 1:numel(pat)
    h1(k) = plot(ax1, pat(k).eps, pat(k).p, ls{k}, ...
        'Color',col(k,:), 'LineWidth',STY.lw);
end
for k = 1:numel(pat)                    % endpoints last, so nothing covers them
    plot(ax1, pat(k).eps(end), pat(k).p(end), 'o', 'MarkerSize',STY.ms, ...
        'MarkerFaceColor',col(k,:), 'MarkerEdgeColor',col(k,:));
    if LABEL_ENDS       % a switch, not dead code
        text(ax1, pat(k).eps(end)+0.025, pat(k).p(end)-0.045, ...
            sprintf('%.3f, %.2f', pat(k).eps(end), pat(k).p(end)), ...
            'Color',col(k,:), 'FontSize',STY.fs-1, 'HorizontalAlignment','left');   %#ok<UNRCH>
    end
end

xlabel(ax1,'Strain \epsilon','FontSize',STY.fsLab);
ylabel(ax1,'Porosity p','FontSize',STY.fsLab);
set(ax1,'XLim',[0 xmax],'YLim',[0 0.8],'XTick',0:0.2:1.2,'YTick',0:0.2:0.8);
lg1 = legend(ax1, h1, short, 'Location','southeast');

%% ---------------------------------- figure 2: normalised rotational energy
[fig2, ax2] = journal_axes(2, STY);

h2 = gobjects(numel(pat),1);
for k = 1:numel(pat)
    h2(k) = plot(ax2, pat(k).eps, pat(k).u, ls{k}, ...
        'Color',col(k,:), 'LineWidth',STY.lw);
end
for k = 1:numel(pat)
    plot(ax2, pat(k).eps(end), pat(k).u(end), 'o', 'MarkerSize',STY.ms, ...
        'MarkerFaceColor',col(k,:), 'MarkerEdgeColor',col(k,:));
    if LABEL_ENDS       % a switch, not dead code
        text(ax2, pat(k).eps(end)+0.025, pat(k).u(end)-0.45, ...
            sprintf('%.3f, %.2f', pat(k).eps(end), pat(k).u(end)), ...
            'Color',col(k,:), 'FontSize',STY.fs-1, 'HorizontalAlignment','left');   %#ok<UNRCH>
    end
end

xlabel(ax2,'Strain \epsilon','FontSize',STY.fsLab);
ylabel(ax2,'Normalised rotational energy u','FontSize',STY.fsLab);
set(ax2,'XLim',[0 xmax],'YLim',[0 8],'XTick',0:0.2:1.2,'YTick',0:2:8);
lg2 = legend(ax2, h2, short, 'Location','northwest');

set([lg1 lg2], 'FontName',STY.font, 'FontSize',STY.fs, 'Box','on', ...
    'EdgeColor',[.2 .2 .2], 'LineWidth',0.5, 'ItemTokenSize',[18 8]);

if SAVE_FIG             % a switch, not dead code
    here = fileparts(mfilename('fullpath'));   %#ok<UNRCH>
    exportgraphics(fig1, fullfile(here,'porosity.pdf'), 'ContentType','vector');
    exportgraphics(fig2, fullfile(here,'rotational_energy.pdf'), 'ContentType','vector');
    fprintf('figures written to %s\n', here);
end

%% ------------------------------------------------------------------ checks
fprintf('fully open state\n');
fprintf('%-20s %8s %8s %8s %10s\n','pattern','theta','strain','porosity','energy');
for k = 1:numel(pat)
    fprintf('%-20s %7.0f%s %8.3f %8.3f %10.3f\n', pat(k).name, ...
        pat(k).theta_max, char(176), pat(k).eps(end), pat(k).p(end), pat(k).u(end));
end

expect = [sqrt(2)-1  0.50;      % squares:   strain, porosity
          1.0        0.75];     % triangles: strain, porosity
got    = [pat(1).eps(end) pat(1).p(end); pat(2).eps(end) pat(2).p(end)];
e_end = max(abs(got(:) - expect(:)));
if e_end < 1e-4, verdict = 'PASS'; else, verdict = 'FAIL'; end
fprintf('\nendpoint check                  max |error| = %.2e  -> %s\n', e_end, verdict);

% lambda here must be the same lambda the deployment code uses
err = 0;
for th_deg = 0:5:90
    S   = forward_deployment('square', 2, 2, 1, th_deg*pi/180);
    err = max(err, abs(S.lambda - sin((pi/2 + th_deg*pi/180)/2)/sin(pi/4)));
end
for th_deg = 0:5:120
    S   = forward_deployment('triangle', 2, 2, 1, th_deg*pi/180);
    err = max(err, abs(S.lambda - sin((pi/3 + th_deg*pi/180)/2)/sin(pi/6)));
end
fprintf('lambda vs forward_deployment    max |error| = %.2e\n', err);

%% ----------------------------------------------------------------- helper
function [fig, ax] = journal_axes(n, STY)
%JOURNAL_AXES  A single-column figure with one styled, gridless axes in it.
fig = figure(n); clf;
set(fig,'Color','w','Units','inches','Position',[1 + 3.7*(n-1), 1, 3.5, 3.0], ...
    'DefaultAxesFontName',STY.font,'DefaultTextFontName',STY.font, ...
    'DefaultAxesFontSize',STY.fs,'DefaultTextFontSize',STY.fs);
ax = axes(fig); hold(ax,'on');
set(ax,'FontName',STY.font,'FontSize',STY.fs,'LineWidth',STY.lwAxis, ...
    'Box','on','TickDir','in','TickLength',[0.018 0.018], ...
    'XColor','k','YColor','k','Layer','top');
grid(ax,'off');
end
