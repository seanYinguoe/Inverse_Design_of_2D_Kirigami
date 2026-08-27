% MAIN  Driver script for the kinematic (stage-1) inverse design of 2D kirigami.
%
%   Reference
%     Ying, Fernando & Dias (2025), "Inverse design of programmable
%     shape-morphing kirigami structures", Int. J. Mech. Sci. 286, 109840.
%     This script implements Section 2.1 "Kinematics" / Fig. 2 of the paper:
%     the constrained kinematic optimisation that produces a cut pattern whose
%     deployed boundary matches a prescribed target shape.
%
%   Pipeline (see README.md for the full workflow)
%     1. fit_initial_guess        - choose the opening angle (and optionally
%                                   the unit length) that starts closest to
%                                   the target, as described in Section 2.1.2.
%     2. tessellation_deployment  - build the m-by-n rotating-squares pattern
%                                   in its compact state for export.
%     3. create_svg_tessellation  - export the compact pattern as a cuttable
%                                   SVG (laser-cutter / COMSOL geometry).
%     4. tessellation_optimization- fmincon over all node coordinates subject
%                                   to the geometric conditions of Section 2.1.
%     5. tessellation_compaction  - fold the optimised deployed pattern back
%                                   to its compact (as-cut) configuration.
%
%   Node numbering used throughout: see create_unit.m.

clear;
clc;
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'function')));

%% Grid
n_rows      = 4;       % number of tessellation units along y
n_cols      = 4;       % number of tessellation units along x
% NOTE both should be EVEN. The symmetry conditions in rigid.m / nonrigid.m
% loop over 1:m/2 and 1:n/2, so an odd grid silently drops the middle row and
% column of symmetry constraints.

%% Target boundary condition and deployability mode
target_shape = 1;      % 1 circle | 2 ellipse | 3 vase | 4 wavy | 5 heart (see shape.m)
shape_size   = 3.2*sqrt(2.0/2.2);   % characteristic size r of the target shape
                                    % (circle radius; ellipse semi-axis a)
deployability = 2;     % 1 rigid-deployable | 2 non-rigid deployable
if deployability == 1; mode = 'rigid'; else; mode = 'nonrigid'; end

%% Initial guess
% The paper picks the opening angle by "minimising the difference between the
% initial configuration and the target shape". Doing that here instead of
% hard-coding an angle is the single cheapest way to make the solve robust.
%
%   FIT_UNIT_LENGTH = false -> unit_length stays fixed, only the angle is fitted
%   FIT_UNIT_LENGTH = true  -> the overall size is fitted too, and the size
%                              constraint is relaxed to an aspect ratio so the
%                              optimiser can keep adjusting it (FreeScale).
FIT_UNIT_LENGTH = true;
unit_length     = 1;   % used as the fixed size, or as the starting size

if FIT_UNIT_LENGTH
    fit_args = {};              % omit the size -> it gets fitted too
else
    fit_args = {unit_length};   % pin the size, fit the angle alone
end
[tessellation_transformed, guess_info] = fit_initial_guess( ...
    n_rows, n_cols, target_shape, shape_size, mode, fit_args{:});
fprintf('initial guess: opening angle %.2f deg, unit length %.4f, boundary RMS %.4f\n', ...
    guess_info.opening_angle*180/pi, guess_info.unit_length, guess_info.residual);

% compact (as-cut) state at the fitted size: opening angle 0, cuts are slits
tessellation_initial = tessellation_deployment(n_rows, n_cols, guess_info.unit_length, 0);

%% Fabrication parameters of the physical cut pattern
ligament_length = 0.30 * 1/n_rows;  % length t of the ligament left at each hinge
cut_width       = 0.50 * ligament_length;   % width w of the cut void

%% Export the compact tessellation to an SVG file. Every cut is shortened by
% ligament_length at its hinge end so the panels stay connected by a small
% ligament; the cut is drawn as a void of width cut_width with filleted ends.
% svg_name = sprintf('tessellation_initial_%dX%d.svg', 2*n_cols, 2*n_rows);
% create_svg_tessellation(tessellation_initial, ligament_length, svg_name, cut_width);

%% Kinematic optimisation
opts = struct();
opts.FreeScale = FIT_UNIT_LENGTH;   % let the optimiser keep tuning the size
% opts.MaxFunctionEvaluations = 400*2*n_rows*n_cols*16;   % default, scales with the grid
% opts.Restarts      = 3;      % warm restarts if fmincon hits a limit
% opts.MinEdgeLength = 0.075;  % floor on panel edges; raise if panels collapse
[tessellation_optimized, info] = tessellation_optimization( ...
    tessellation_transformed, target_shape, shape_size, deployability, opts);

if ~info.converged
    warning('main:notConverged', ...
        ['The optimisation did not converge (max|ceq| = %.2e). The plots below ' ...
         'will look like a tessellation but do NOT satisfy the geometric ' ...
         'conditions, and the compacted state will not be a proper rectangle. ' ...
         'See README section "Making the optimiser converge".'], info.maxceq);
end

% fold the optimised deployed pattern back to its compact, as-cut state
tessellation_compacted = tessellation_compaction(tessellation_optimized);

%% Plot the results
figure(2); clf;
subplot(2,2,1);
title('intial kirigami tessellation');
plot_tessellation(tessellation_initial);
axis off

subplot(2,2,2);
title('deployed kirigami tessellation');
plot_tessellation(tessellation_transformed);
plot_boundary(target_shape, shape_size);
axis off

subplot(2,2,3);
title(sprintf('optimized (max|ceq| = %.1e)', info.maxceq));
plot_tessellation(tessellation_optimized);
plot_boundary(target_shape, shape_size);
axis off

subplot(2,2,4);
title('valid compacted kirigami tessellation');
plot_tessellation(tessellation_compacted);
axis off
