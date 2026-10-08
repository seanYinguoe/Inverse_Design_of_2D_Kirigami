% FORWARD_ANISOTROPY  Anisotropic rotating-units kirigami: rectangles and parallelograms.
%
%   SCRIPT. The anisotropic generalisation of FORWARD_ANALYSIS. Replacing the
%   square panels of the regular mechanism by parallelograms of edge lengths
%   a, b and included angle gamma turns the single isotropic stretch into two
%   different principal stretches, while keeping the one-degree-of-freedom
%   character of the motion. Rotating squares (a = b, gamma = 90 deg) is
%   recovered exactly and is the ONLY isotropic member of the family.
%
%   MODEL  (derived in PARALLELOGRAM_METRICS)
%     With u, v the panel edge vectors and alpha the PANEL rotation,
%
%       s(alpha) = u cos(alpha) - Jv sin(alpha)      deployed lattice vectors
%       t(alpha) = v cos(alpha) + Ju sin(alpha)
%       F(alpha) = cos(alpha) I - sin(alpha) M,      M = J U J U^(-1)
%
%     M is symmetric, so F is symmetric and its principal frame is FIXED along
%     the whole deployment path. det M = 1 and tr M = -(a^2+b^2)/(ab sin gamma),
%     hence the two anisotropy factors satisfy
%
%       mu1 * mu2 = 1        mu1 + mu2 = (a^2 + b^2)/(a b sin gamma)
%
%     and the principal stretches, porosity and Poisson's ratio are
%
%       lambda_i = cos(alpha) + mu_i sin(alpha)
%       p        = 1 - 1/(lambda1 lambda2)
%       nu       = -dln(lambda2)/dln(lambda1)        on the principal axes
%
%     Edge ratio and panel angle feed the SAME scalar mu, so different (a/b,
%     gamma) pairs can be kinematically identical.
%
%   FIGURES  (all in the FORWARD_ANALYSIS style)
%     1  unit kinematics - compact and deployed sheet for two design points,
%        with the deployed lattice vectors and one void marked
%     2  principal stretches lambda_1, lambda_2 against alpha
%     3  design map of mu over (aspect ratio, panel angle)
%
%   See also PARALLELOGRAM_METRICS, PARALLELOGRAM_DEPLOYMENT, FORWARD_ANALYSIS,
%   PLOT_FORWARD, JOURNAL_AXES.

clear;
clc;
addpath(fileparts(fileparts(fileparts(mfilename('fullpath'))))); setup_project;

%% =========================================================== CONFIG BLOCK
% Every geometric constant in this script lives here and nowhere else.
CFG.shape = 'triangle';        % 'parallelogram' | 'triangle'

CFG.b     = 1;                 % reference panel edge length
CFG.grid  = [3 3];             % panels drawn in the kinematics figure [m n]

switch lower(CFG.shape)
    case 'parallelogram'       % four panels meet at a vertex; isotropic = square
        CFG.ratios    = [1.5 1.5];      % a/b at the two design points
        CFG.gammas    = [90  65];       % gamma [deg] at the same two
        CFG.iso       = [1 90];         % the family's single isotropic point
        CFG.iso_name  = 'rotating squares';
        CFG.map_gamma = [35 90];
        CFG.levels    = [1.2 1.5 2 2.5 3 4];
        CFG.ylim      = [1 2.1];
    case 'triangle'            % six panels meet; isotropic = equilateral
        CFG.ratios    = [1.3 1.3];
        CFG.gammas    = [60  45];
        CFG.iso       = [1 60];
        CFG.iso_name  = 'rotating equilateral triangles';
        CFG.map_gamma = [30 90];
        CFG.levels    = [2 2.5 3 4 5 6];
        CFG.ylim      = [1 3.1];
    otherwise
        error('forward_anisotropy:shape','CFG.shape must be parallelogram or triangle');
end
CFG.map_ratio  = [1 3];         % design-map range of a/b
CFG.show_frac  = 2/3;           % where along the sweep the deployed pictures sit.
                                % Set as a FRACTION, not a fixed angle: the
                                % deployed cell angle is degenerate at exactly
                                % alpha = 45 deg for parallelograms and 30 deg
                                % for triangles (s.t = 0 there for EVERY panel
                                % of the family), so a fixed angle can land on
                                % the one rotation where the annotation carries
                                % no information.

% End of the sweep: the rotation at which the family's ISOTROPIC member is
% fully open, atan(sqrt(k)) - 45 deg for parallelograms, 60 deg for triangles.
% A sweep convention, not an admissibility limit; see ROTATING_UNIT_METRICS.
REF = rotating_unit_metrics(CFG.shape, CFG.b, CFG.b, CFG.gammas(1)*pi/180);
CFG.alpha_max = REF.alpha_iso*180/pi;
CFG.k         = REF.k;
CFG.alpha_show = round(CFG.show_frac*CFG.alpha_max);
CFG.labels    = arrayfun(@(g,r) sprintf('\\gamma = %g\\circ,  a/b = %g', g, r), ...
                         CFG.gammas, CFG.ratios, 'UniformOutput', false);

STY = struct('font','Helvetica', 'fs',9, 'fsLab',10, 'lw',1.6, ...
             'lwAxis',0.9, 'ms',5);

SAVE_FIG = false;                   % write vector PDFs next to this script
OUTDIR   = fileparts(mfilename('fullpath'));

C   = forward_colors();
col = C([1 2],:);                   % blue = point 1, orange = point 2
ls  = {'-','--'};

%% ================================================== figure 1: kinematics
[fig1, AX] = journal_axes(1, STY, [5.6 4.6], [2 2]);
state = {'compact','deployed'};

for r = 1:2
    a  = CFG.ratios(r)*CFG.b;
    g  = CFG.gammas(r)*pi/180;
    for cc = 1:2
        al = (cc-1)*CFG.alpha_show*pi/180;
        S  = rotating_unit_deployment(CFG.shape, CFG.grid(1), CFG.grid(2), a, CFG.b, g, al);
        ax = AX(r,cc);
        plot_forward(S, 'Parent', ax);

        if cc == 2
            % Lattice vectors s and t at true length, anchored at the centre of
            % the middle panel - s and t connect panel centres, so the arrows
            % land on the neighbouring panels. The void boundary is already
            % drawn: it is made of panel edges, which PLOT_FORWARD colours.
            o = mean(S.panels{ceil(numel(S.panels)/2)}, 1);
            quiver(ax, [o(1) o(1)], [o(2) o(2)], [S.s(1) S.t(1)], [S.s(2) S.t(2)], 0, ...
                'Color',[.1 .1 .1], 'LineWidth',1.3, 'MaxHeadSize',0.35);
            text(ax, o(1), o(2)-0.14*S.scale, ...
                sprintf('%.1f\\circ', abs(S.cellangle)*180/pi), ...
                'Color',[.1 .1 .1], 'FontSize',STY.fs-1, ...
                'HorizontalAlignment','center','VerticalAlignment','top');
        end
        title(ax, sprintf('%s  (\\alpha = %.0f\\circ)', state{cc}, al*180/pi), ...
            'FontSize',STY.fs, 'FontWeight','normal');
    end
    % PLOT_FORWARD switches the axes off, so ylabel would never appear;
    % place the row label as rotated text in normalised coordinates instead.
    text(AX(r,1), -0.06, 0.5, CFG.labels{r}, 'Units','normalized', ...
        'Rotation',90, 'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', 'FontSize',STY.fsLab);
end

%% ========================================= figure 2: principal stretches
[fig2, ax2] = journal_axes(2, STY);

al  = linspace(0, CFG.alpha_max, 401)*pi/180;
h2  = gobjects(2,1);
for r = 1:2
    M = rotating_unit_metrics(CFG.shape, CFG.ratios(r)*CFG.b, CFG.b, CFG.gammas(r)*pi/180, al);
    h2(r) = plot(ax2, al*180/pi, M.lambda(1,:), ls{r}, ...
        'Color',col(1,:), 'LineWidth',STY.lw);
    plot(ax2, al*180/pi, M.lambda(2,:), ls{r}, ...
        'Color',col(2,:), 'LineWidth',STY.lw);
    fprintf('design point %d: a/b = %.2f, gamma = %2d deg  ->  mu = %.4f\n', ...
        r, CFG.ratios(r), CFG.gammas(r), M.mu(1));
end
xlabel(ax2, 'Panel rotation \alpha (deg)', 'FontSize',STY.fsLab);
ylabel(ax2, 'Principal stretches \lambda_{1,2}', 'FontSize',STY.fsLab);
set(ax2, 'XLim',[0 CFG.alpha_max], 'YLim',CFG.ylim, ...
         'XTick',0:10:CFG.alpha_max);
lg2 = legend(ax2, h2, ...
    arrayfun(@(r) sprintf('\\gamma = %g\\circ,  \\mu = %.2f', CFG.gammas(r), ...
        rotating_unit_metrics(CFG.shape,CFG.ratios(r)*CFG.b,CFG.b,CFG.gammas(r)*pi/180).mu(1)), ...
        1:2, 'UniformOutput',false), ...
    'Location','northwest');
M1 = rotating_unit_metrics(CFG.shape, CFG.ratios(1)*CFG.b, CFG.b, ...
                           CFG.gammas(1)*pi/180, CFG.alpha_max*pi/180);
text(ax2, CFG.alpha_max-1.5, M1.lambda(1)+0.07, '\lambda_1', ...
    'Color',col(1,:), 'FontSize',STY.fsLab, 'HorizontalAlignment','right');
text(ax2, CFG.alpha_max-1.5, M1.lambda(2)+0.07, '\lambda_2', ...
    'Color',col(2,:), 'FontSize',STY.fsLab, 'HorizontalAlignment','right');

%% =============================================== figure 3: mu design map
[fig3, ax3] = journal_axes(3, STY);

rr = linspace(CFG.map_ratio(1), CFG.map_ratio(2), 260);
gg = linspace(CFG.map_gamma(1), CFG.map_gamma(2), 260);
[RR, GG] = meshgrid(rr, gg);
% H = sum(edge^2)/(2 * 2*area), vectorised over the grid; mu = H + sqrt(H^2-k).
% Spot-checked against ROTATING_UNIT_METRICS in the checks block below.
switch lower(CFG.shape)
    case 'parallelogram'
        HH = (RR.^2 + 1) ./ (2*RR.*sind(GG));
    case 'triangle'
        HH = (RR.^2 + 1 - RR.*cosd(GG)) ./ (RR.*sind(GG));
end
MU = HH + sqrt(max(HH.^2 - CFG.k, 0));

contourf(ax3, RR, GG, MU, 24, 'LineStyle','none');
[Cn, hn] = contour(ax3, RR, GG, MU, CFG.levels, ...
    'LineColor',[.25 .25 .25], 'LineWidth',0.6);
clabel(Cn, hn, 'FontSize',STY.fs-2, 'Color',[.25 .25 .25], ...
    'LabelSpacing',260, 'FontName',STY.font);
colormap(ax3, parula);
cb = colorbar(ax3);
cb.Label.String = '\mu';
set(cb,'FontName',STY.font,'FontSize',STY.fs,'LineWidth',STY.lwAxis);

plot(ax3, CFG.iso(1), CFG.iso(2), 's', 'MarkerFaceColor',col(2,:), ...
    'MarkerEdgeColor','w', 'MarkerSize',STY.ms+2, 'LineWidth',0.8);
text(ax3, CFG.iso(1)+0.10, CFG.iso(2)-3, {CFG.iso_name,'(only isotropic point)'}, ...
    'FontSize',STY.fs-1, 'Color','w', 'VerticalAlignment','top');
for r = 1:2
    plot(ax3, CFG.ratios(r), CFG.gammas(r), 'o', 'MarkerFaceColor',col(2,:), ...
        'MarkerEdgeColor','w', 'MarkerSize',STY.ms, 'LineWidth',0.8);
end

xlabel(ax3, 'Aspect ratio a/b', 'FontSize',STY.fsLab);
ylabel(ax3, 'Panel angle \gamma (deg)', 'FontSize',STY.fsLab);
set(ax3, 'XLim',CFG.map_ratio, 'YLim',CFG.map_gamma, 'Layer','top');

set(lg2, 'FontName',STY.font, 'FontSize',STY.fs, 'Box','on', ...
    'EdgeColor',[.2 .2 .2], 'LineWidth',0.5, 'ItemTokenSize',[18 8]);

%% ------------------------------------------------------------------ export
if SAVE_FIG             % a switch, not dead code
    tag = CFG.shape;   %#ok<UNRCH>
    exportgraphics(fig1, fullfile(OUTDIR,[tag '_unit_kinematics.pdf']), 'ContentType','vector');
    exportgraphics(fig2, fullfile(OUTDIR,[tag '_principal_stretches.pdf']), 'ContentType','vector');
    exportgraphics(fig3, fullfile(OUTDIR,[tag '_design_map.pdf']), 'ContentType','vector');
    fprintf('figures written to %s\n', OUTDIR);
end

%% ------------------------------------------------------------------ checks
fprintf('\n%-40s %s\n','check','max |error|');

% (1) regular-case regression: the isotropic member of this family must
%     reproduce the existing regular-case code exactly.
switch lower(CFG.shape)
    case 'parallelogram'
        ref = @(th) square_deployment(2,2,1,th);        % existing code
        gi  = pi/2;   thmax = 90;
    case 'triangle'
        ref = @(th) triangle_deployment(2,2,1,th);      % existing code
        gi  = pi/3;   thmax = 120;
end
eL = 0; eN = 0; eP = 0; eG = 0;
for th = (0:5:thmax)*pi/180
    Rg = ref(th);
    Mm = rotating_unit_metrics(CFG.shape, 1, 1, gi, th/2);
    eL = max(eL, max(abs(Mm.lambda - Rg.lambda)));
    eP = max(eP, abs(Mm.p - (1 - Rg.lambda^-2)));
    if th > 0, eN = max(eN, abs(Mm.nu + 1)); end
    eG = max(eG, rotating_unit_deployment(CFG.shape,3,3,1,1,gi,th/2).hinge_gap);
end
fprintf('%-40s %.2e\n','isotropic limit: lambda_i',  eL);
fprintf('%-40s %.2e\n','isotropic limit: nu = -1',   eN);
fprintf('%-40s %.2e\n','isotropic limit: porosity',  eP);
fprintf('%-40s %.2e\n','isotropic limit: hinge gap', eG);

% (2) the vectorised design-map formula must agree with the metrics function
eM = 0;
for r = linspace(CFG.map_ratio(1)+0.01, CFG.map_ratio(2), 7)
    for g = linspace(CFG.map_gamma(1), CFG.map_gamma(2), 7)
        Mm = rotating_unit_metrics(CFG.shape, r*CFG.b, CFG.b, g*pi/180);
        switch lower(CFG.shape)
            case 'parallelogram', Hc = (r^2+1)/(2*r*sind(g));
            case 'triangle',      Hc = (r^2+1-r*cosd(g))/(r*sind(g));
        end
        eM = max(eM, abs(Mm.mu(1) - (Hc + sqrt(max(Hc^2-CFG.k,0)))));
    end
end
fprintf('%-40s %.2e\n','design map vs rotating_unit_metrics', eM);

% (3) geometry of the drawn design points, and the exact overlap test
fprintf('\n%-14s %8s %11s %11s %13s\n', ...
    'design point','mu','angle(s,t)','180-dir(t)','overlap area');
for r = 1:2
    Sd = rotating_unit_deployment(CFG.shape, CFG.grid(1), CFG.grid(2), ...
        CFG.ratios(r)*CFG.b, CFG.b, CFG.gammas(r)*pi/180, CFG.alpha_show*pi/180);
    fprintf('%-14d %8.4f %10.2f%s %10.2f%s %13.2e\n', r, Sd.mu(1), ...
        abs(Sd.cellangle)*180/pi, char(176), ...
        180 - atan2d(Sd.t(2),Sd.t(1)), char(176), Sd.overlap);
end

% (4) nu changes sign at alpha = atan(mu2): the minor direction starts to
%     contract while the major one still expands, so the sheet stops being
%     auxetic partway along the path.
fprintf('\n');
for r = 1:2
    Mm = rotating_unit_metrics(CFG.shape, CFG.ratios(r)*CFG.b, CFG.b, ...
                               CFG.gammas(r)*pi/180, al);
    fprintf(['design point %d: mu = %.4f, nu(0) = %+.4f, ' ...
             'nu changes sign at alpha = %.2f deg\n'], ...
        r, Mm.mu(1), Mm.nu(1), atan(Mm.mu(2))*180/pi);
end
