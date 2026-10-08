function [tessellation_guess, best] = fit_initial_guess(m, n, s, r, mode, unit_length, angle_range)
%FIT_INITIAL_GUESS  Pick the opening angle (and size) that starts closest to the target.
%
%   [tessellation_guess, best] = fit_initial_guess(m, n, s, r)
%   [tessellation_guess, best] = fit_initial_guess(m, n, s, r, mode)
%   [tessellation_guess, best] = fit_initial_guess(m, n, s, r, mode, unit_length)
%
%   Implements the step described in Section 2.1.2 of the paper: "to reduce
%   the optimisation time, we minimise the difference between the initial
%   configuration and the target shape by controlling the opening angle of
%   the square tessellation". The code previously used a hard-coded opening
%   angle, which left the optimiser starting far from a feasible point.
%
%   The search is over the two parameters of a UNIFORM deployment:
%     xi           - the opening angle, the single degree of freedom of the
%                    rotating-squares mechanism
%     unit_length  - the size of one unit, which sets the overall scale of
%                    the sheet and therefore how well it can fill the target
%   Both are cheap to evaluate (no constraints, just tessellation_deployment),
%   so a coarse grid followed by a Nelder-Mead refinement is plenty.
%
%   INPUTS
%     m, n        : number of units along y and x
%     s           : target shape selector, see shape.m
%     r           : characteristic size of the target shape
%     mode        : 'rigid' (default) or 'nonrigid' - passed to
%                   boundary_residual, only matters for shape 3
%     unit_length : fix the unit length to this value and search the opening
%                   angle alone. Omit or leave empty to optimise both.
%     angle_range : [lo hi] bounds on the opening angle, default
%                   [0.10  pi/2-0.10] i.e. about 5.7 to 84.3 degrees.
%                   The bounds deliberately exclude the two ENDS of the
%                   mechanism's travel. Minimising boundary mismatch alone
%                   often drives xi to 0 (fully closed) or pi/2, but those
%                   are exactly the configurations where every void has zero
%                   area, so the non-overlap inequalities sit at equality.
%                   An interior-point solver started on its own constraint
%                   boundary has nowhere to step, so a slightly worse start
%                   that is comfortably feasible converges far better.
%
%   OUTPUTS
%     tessellation_guess : m-by-n cell array, the best uniform deployment
%     best : struct with fields
%              opening_angle - the fitted xi in radians
%              unit_length   - the fitted unit length
%              residual      - RMS boundary mismatch of the guess
%              searched      - 'angle' or 'angle+length'
%
%   Feed tessellation_guess straight to tessellation_optimization.

if nargin < 5 || isempty(mode);        mode = 'rigid'; end
if nargin < 6;                         unit_length = []; end
if nargin < 7 || isempty(angle_range); angle_range = [0.10, pi/2 - 0.10]; end

fixL = ~isempty(unit_length);
if fixL
    L0 = unit_length;
    best.searched = 'angle';
else
    L0 = 1;
    best.searched = 'angle+length';
end

% --- objective: RMS of the boundary residuals of a uniform deployment ----
    function v = cost(xi, L)
        % keep the search inside the physically meaningful range
        if ~isfinite(xi) || ~isfinite(L) || L <= 1e-3 || ...
                xi < angle_range(1) || xi > angle_range(2)
            v = 1e6; return
        end
        try
            T = tessellation_deployment(m, n, L, xi);
            [~, ceq] = boundary_residual(T, s, r, m, n, mode);
        catch
            v = 1e6; return
        end
        if isempty(ceq) || any(~isfinite(ceq))
            v = 1e6; return
        end
        v = sqrt(mean(ceq.^2));
    end

% --- coarse grid, so the refinement cannot start in a bad basin ----------
xiGrid = linspace(angle_range(1), angle_range(2), 60);
if fixL
    LGrid = L0;
else
    LGrid = L0 * linspace(0.4, 2.5, 40);
end
bestV = inf; bestXi = xiGrid(1); bestL = LGrid(1);
for a = 1:numel(xiGrid)
    for b = 1:numel(LGrid)
        v = cost(xiGrid(a), LGrid(b));
        if v < bestV
            bestV = v; bestXi = xiGrid(a); bestL = LGrid(b);
        end
    end
end

% --- local refinement ----------------------------------------------------
opts = optimset('Display','off','TolX',1e-8,'TolFun',1e-10,'MaxFunEvals',2000);
if fixL
    f = @(z) cost(z, L0);
    [z, v] = fminsearch(f, bestXi, opts);
    bestXi = z; bestL = L0; bestV = v;
else
    f = @(z) cost(z(1), z(2));
    [z, v] = fminsearch(f, [bestXi bestL], opts);
    if v < bestV
        bestXi = z(1); bestL = z(2); bestV = v;
    end
end

tessellation_guess = tessellation_deployment(m, n, bestL, bestXi);
best.opening_angle = bestXi;
best.unit_length   = bestL;
best.residual      = bestV;
best.at_bound      = abs(bestXi - angle_range(1)) < 1e-3 || ...
                     abs(bestXi - angle_range(2)) < 1e-3;
if best.at_bound
    warning('fit_initial_guess:atBound', ...
        ['The best opening angle sits on the search bound (%.1f deg). The ' ...
         'target may be a poor match for a %d-by-%d uniform grid; consider a ' ...
         'different grid size, or widen angle_range if you accept a start ' ...
         'close to the closed or fully open state.'], bestXi*180/pi, m, n);
end
end
