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
%     1. tessellation_deployment  - build the m-by-n rotating-squares pattern,
%                                   either compact (opening angle 0) or partly
%                                   deployed (opening angle > 0).
%     2. create_svg_tessellation  - export the compact pattern as a cuttable
%                                   SVG (laser-cutter / COMSOL geometry).
%     3. tessellation_optimization- fmincon over all node coordinates subject
%                                   to the geometric conditions of Section 2.1.
%     4. tessellation_compaction  - fold the optimised deployed pattern back
%                                   to its compact (as-cut) configuration.
%
%   Node numbering used throughout: see create_unit.m.

clear;
clc;
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'function')));

%% Initial configuration of the kirigami tessellation
n_rows       = 4;      % number of tessellation units along y
n_cols       = 4;      % number of tessellation units along x
unit_length  = 1;      % side length of one unit (= 2 quadrant squares)
opening_angle = pi/6;  % rotation angle of the squares in the deployed state

% compact (as-cut) state: opening angle 0, all cuts are zero-width slits
tessellation_initial     = tessellation_deployment(n_rows, n_cols, unit_length, 0);
% uniformly deployed state: used as the initial guess for the optimiser
tessellation_transformed = tessellation_deployment(n_rows, n_cols, unit_length, opening_angle);

%% Target boundary condition and deployability mode
target_shape = 1;      % 1 circle | 2 ellipse | 3 vase | 4 wavy | 5 heart (see shape.m)
shape_size   = 3.2*sqrt(2.0/2.2);   % characteristic size r of the target shape
                                    % (circle radius; ellipse semi-axis a)
deployability = 2;     % 1 rigid-deployable | 2 non-rigid deployable

%% Fabrication parameters of the physical cut pattern
ligament_length = 0.30 * 1/n_rows;  % length t of the ligament left at each hinge
cut_width       = 0.50 * ligament_length;   % width w of the cut void

%% Export the compact tessellation to an SVG file. Every cut is shortened by
% ligament_length at its hinge end so the panels stay connected by a small
% ligament; the cut is drawn as a void of width cut_width with filleted ends.
% svg_name = sprintf('tessellation_initial_%dX%d.svg', 2*n_cols, 2*n_rows);
% create_svg_tessellation(tessellation_initial, ligament_length, svg_name, cut_width);

%% Kinematic optimisation (uncomment to run - takes minutes for large grids)
tessellation_optimized = tessellation_optimization(tessellation_transformed, ...
                                                    target_shape, shape_size, deployability);
% % fold the optimised deployed pattern back to its compact, as-cut state
tessellation_compacted = tessellation_compaction(tessellation_optimized);

%% Plot the results
figure(2);
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
title('optimized kirigami tessellation');
plot_tessellation(tessellation_optimized);
plot_boundary(target_shape, shape_size);
axis off

subplot(2,2,4);
title('valid compacted kirigami tessellation');
plot_tessellation(tessellation_compacted);
axis off
