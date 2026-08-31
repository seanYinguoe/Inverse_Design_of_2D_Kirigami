function S = forward_deployment(unit_type, m, n, unit_length, theta)
%FORWARD_DEPLOYMENT  Forward kinematics of a uniform kirigami tessellation.
%
%   S = forward_deployment(unit_type, m, n, unit_length, theta)
%
%   Opposite of the inverse design in MAIN.M: instead of solving for a cut
%   pattern that reaches a prescribed target shape, this simply sweeps the
%   single degree of freedom of a periodic rotating-units mechanism and
%   reports the resulting configuration.
%
%   INPUTS
%     unit_type   : 'square'   - rotating squares,   theta in [0, pi/2]
%                   'triangle' - rotating triangles, theta in [0, 2*pi/3]
%     m, n        : number of cells along y and x. A square cell holds 4
%                   quadrant panels, a triangular cell holds 2, so the patch
%                   has 4*m*n or 2*m*n rigid panels.
%     unit_length : side of one square unit / side of one triangle
%     theta       : HINGE opening angle in radians, measured between the two
%                   panels that meet at a hinge. 0 is the compact, as-cut
%                   state; theta_max is the fully open one. Each panel itself
%                   turns by theta/2.
%
%   OUTPUT  S, a struct describing one configuration:
%     .type       'square' | 'triangle'
%     .theta      the opening angle it was evaluated at [rad]
%     .theta_max  the fully open angle for this unit type [rad]
%     .panels     P-by-1 cell, each K-by-2 list of panel corners (K = 4 or 3)
%     .edges.seg  E-by-4 [x1 y1 x2 y2], one row per shared (cut) edge. A cut
%                 contributes TWO rows, one per panel, so you can watch the
%                 pair separate as theta grows.
%     .edges.cid  E-by-1 colour index into FORWARD_COLORS, constant in theta
%     .hinges     H-by-2 hinge coordinates
%     .hinge_gap  largest distance between two corners that should coincide
%                 at a hinge - a numerical check that the motion is rigid,
%                 should stay at round-off
%     .ncolors    4 (square) or 6 (triangle)
%     .scale      side of one rigid panel, for line widths / marker sizes
%     .lambda     linear expansion ratio with respect to the compact state
%     .view_span  width of the fully open patch, used to weight line widths
%
%   The patch is centred on the origin at every theta, so a sweep animates in
%   place. See FORWARD_MAIN for a driver, PLOT_FORWARD for the drawing.
%
%   See also FORWARD_MAIN, PLOT_FORWARD, SQUARE_DEPLOYMENT,
%   TRIANGLE_DEPLOYMENT, TESSELLATION_DEPLOYMENT.

switch lower(unit_type)
    case {'square','squares'}
        S = square_deployment(m, n, unit_length, theta);
    case {'triangle','triangles','triangular'}
        S = triangle_deployment(m, n, unit_length, theta);
    otherwise
        error('forward_deployment:unitType', ...
            'unit_type must be ''square'' or ''triangle'', got ''%s''.', unit_type);
end

% centre the patch on the origin so that a theta sweep animates in place
allv = vertcat(S.panels{:});
c    = (min(allv,[],1) + max(allv,[],1))/2;
for p = 1:numel(S.panels)
    S.panels{p} = S.panels{p} - c;
end
S.edges.seg = S.edges.seg - [c c];
S.hinges    = S.hinges - c;
end
