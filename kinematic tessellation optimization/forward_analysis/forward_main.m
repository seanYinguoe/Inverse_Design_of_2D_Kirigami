% FORWARD_MAIN  Forward (kinematic) deployment of a uniform kirigami tessellation.
%
%   SCRIPT. The forward counterpart of MAIN.M: MAIN.M solves the inverse
%   problem (find the cut pattern whose deployed boundary matches a target
%   shape), while this script just turns the single degree of freedom of a
%   periodic rotating-units mechanism and shows what comes out.
%
%   Set the four knobs below and run:
%
%     unit_type     'square'   - rotating squares,   0 <= theta <= 90 deg
%                   'triangle' - rotating triangles, 0 <= theta <= 120 deg
%     n_rows/n_cols size of the patch
%     unit_length   size of one rigid panel
%     opening_angle the degree of freedom, in DEGREES
%
%   WHAT IS DRAWN
%     * the rigid panels, as pale filled polygons;
%     * every shared edge - the cuts, i.e. the edges that lie on top of each
%       other in the compact state and peel apart as the sheet opens - in a
%       colour that it keeps for the whole motion. Four colours are enough
%       for the square pattern and six for the triangular one, because that
%       is how many cuts meet at one rotation centre: every void of the
%       deployed sheet is bounded by exactly one edge of each colour, so you
%       can read straight off the picture where every edge of the compact
%       pattern ended up;
%     * every hinge as a red dot. Hinges are the only points where two panels
%       stay connected; in the compact state several of them sit on top of
%       each other and they separate as the sheet opens.
%
%   See also FORWARD_DEPLOYMENT, PLOT_FORWARD, MAIN.

clear;
clc;
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'function')));

%% ------------------------------------------------------------------ knobs
unit_type     = 'triangle';   % 'square' | 'triangle'
n_rows        = 3;            % cells along y
n_cols        = 3;            % cells along x
unit_length   = 1;            % side of one square unit / of one triangle

opening_angle = 120;           % DEGREES. 0 = compact (as cut)
                              %   square   : 0 .. 90   (90  = fully open)
                              %   triangle : 0 .. 120  (120 = fully open)

ANIMATE  = false;             % play the whole 0 -> fully open sweep
N_FRAMES = 45;                % frames in that sweep
SAVE_GIF = false;             % write the sweep to <unit_type>_deployment.gif

%% ------------------------------------------------------- the two snapshots
theta   = opening_angle*pi/180;
S_open  = forward_deployment(unit_type, n_rows, n_cols, unit_length, theta);
S_flat  = forward_deployment(unit_type, n_rows, n_cols, unit_length, 0);
S_full  = forward_deployment(unit_type, n_rows, n_cols, unit_length, S_open.theta_max);

fprintf('%s tessellation, %d x %d cells, %d rigid panels\n', ...
    S_open.type, n_rows, n_cols, numel(S_open.panels));
fprintf('opening angle   : %.1f deg  (fully open at %.0f deg)\n', ...
    opening_angle, S_open.theta_max*180/pi);
fprintf('panel rotation  : %.1f deg\n', opening_angle/2);
fprintf('expansion       : %.4f in length, %.4f in area\n', ...
    S_open.lambda, S_open.lambda^2);
fprintf('cuts / hinges   : %d / %d\n', ...
    size(S_open.edges.seg,1)/2, size(S_open.hinges,1));
fprintf('hinge closure   : %.2e  (should be round-off - the motion is rigid)\n', ...
    S_open.hinge_gap);

% common axis limits, taken from the fully open state so nothing jumps
allv = vertcat(S_full.panels{:});
lo   = min(allv,[],1);   hi = max(allv,[],1);
pad  = 0.05*max(hi-lo);
lim  = [lo(1)-pad hi(1)+pad lo(2)-pad hi(2)+pad];

%% ------------------------------------------------------------------- plot
figure(1); clf;

ax1 = subplot(1,2,1);
plot_forward(S_flat, 'Parent', ax1);
axis(ax1, lim);
title(ax1, sprintf('compact  (\\theta = 0\\circ)'));

ax2 = subplot(1,2,2);
plot_forward(S_open, 'Parent', ax2);
axis(ax2, lim);
title(ax2, sprintf('deployed  (\\theta = %.0f\\circ,  %.2f\\times in length)', ...
    opening_angle, S_open.lambda));

sgtitle(sprintf('forward deployment of a %s tessellation, %d \\times %d cells', ...
    S_open.type, n_rows, n_cols));

%% -------------------------------------------------------------- animation
if ANIMATE
    figure(2); clf;   %#ok<UNRCH> ANIMATE is a switch, not dead code
    ax = axes;
    thetas = linspace(0, S_open.theta_max, N_FRAMES);
    gif_name = sprintf('%s_deployment.gif', S_open.type);
    for f = 1:numel(thetas)
        Sf = forward_deployment(unit_type, n_rows, n_cols, unit_length, thetas(f));
        cla(ax);
        plot_forward(Sf, 'Parent', ax);
        axis(ax, lim);
        title(ax, sprintf('\\theta = %5.1f\\circ', thetas(f)*180/pi));
        drawnow;
        if SAVE_GIF
            frame = getframe(gcf);
            [idx, cm] = rgb2ind(frame2im(frame), 256);
            if f == 1
                imwrite(idx, cm, gif_name, 'gif', 'LoopCount', inf, 'DelayTime', 0.05);
            else
                imwrite(idx, cm, gif_name, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);
            end
        end
    end
    if SAVE_GIF
        fprintf('animation written to %s\n', gif_name);
    end
end
