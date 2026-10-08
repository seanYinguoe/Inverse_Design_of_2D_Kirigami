function [c, ceq] = boundary_residual(tessellation, s, r, m, n, mode)
%BOUNDARY_RESIDUAL  Boundary-shape conditions of the deployed tessellation.
%
%   [c, ceq] = boundary_residual(tessellation, s, r, m, n, mode)
%
%   Single source of truth for Eq. (7) of the paper, "boundary shape
%   condition". rigid.m and nonrigid.m call this for their boundary block,
%   and fit_initial_guess.m calls it to choose the opening angle and unit
%   length that start the optimiser closest to the target.
%
%   INPUTS
%     tessellation : m-by-n cell array of 16-by-2 node lists
%     s            : target shape. Either a built-in code 1..5 (see shape.m),
%                    or a target struct from make_target describing ANY
%                    closed curve - an implicit function or a point list.
%                    See make_target for how to build one.
%     r            : characteristic size of the target shape
%     m, n         : number of units along y and x
%     mode         : 'rigid' or 'nonrigid'. Only shape 3 (vase) differs -
%                    the two constraint files historically wrote that one
%                    differently, and both forms are preserved here.
%
%   OUTPUTS
%     c   : inequality residuals, must be <= 0 (nodes inside the shape)
%     ceq : equality residuals, must be 0 (boundary nodes on the shape)
%
%   BOUNDARY NODES, per unit (see create_unit.m)
%     left edge   14, 9      right edge  3, 8
%     bottom edge 12, 5      top edge    15, 2
%
%   For shapes 2-4 the left and right edges are the loading grips: they are
%   clamped straight at x = +-HALF_WIDTH instead of following the curve.
%
%   The residual ORDER is part of the contract - rigid.m and nonrigid.m
%   append these to their running c / ceq vectors, so changing the order
%   changes nothing numerically but do keep it stable for readability.

if nargin < 6 || isempty(mode)
    mode = 'rigid';
end

% --- general target: any closed curve ------------------------------------
% Built-in codes keep their original hard-coded treatment below, so published
% results reproduce exactly. Anything built by make_target goes through the
% general path, which is simply Eq. (7) applied literally: every boundary node
% must sit on the curve, every other node must stay inside it.
if isstruct(s)
    [c, ceq] = general_target(tessellation, s, m, n);
    return
end

HALF_WIDTH = 2.5;   % x position of the clamped edges; see plot_boundary.m

c   = [];
ceq = [];

switch s

    case 1   % ---- circle -------------------------------------------------
        % every node must lie inside the circle
        c = zeros(1, m*n*16);
        q = 1;
        for i = 1:m
            for j = 1:n
                for l = 1:16
                    c(q) = length_calculate([tessellation{i,j}(l,:);[0,0]]) - r;
                    q = q + 1;
                end
            end
        end
        ceq = zeros(1, 4*m + 4*n);
        k = 1;
        for i = 1:m      % left and right boundary nodes
            ceq(k)   = shape(s, tessellation{i,1}(14,:), r);
            ceq(k+1) = shape(s, tessellation{i,1}(9,:),  r);
            ceq(k+2) = shape(s, tessellation{i,n}(3,:),  r);
            ceq(k+3) = shape(s, tessellation{i,n}(8,:),  r);
            k = k + 4;
        end
        for j = 1:n      % bottom and top boundary nodes
            ceq(k)   = shape(s, tessellation{1,j}(12,:), r);
            ceq(k+1) = shape(s, tessellation{1,j}(5,:),  r);
            ceq(k+2) = shape(s, tessellation{m,j}(15,:), r);
            ceq(k+3) = shape(s, tessellation{m,j}(2,:),  r);
            k = k + 4;
        end

    case 2   % ---- ellipse ------------------------------------------------
        c = zeros(1, 2*m*n*16);
        q = 1;
        for i = 1:m
            for j = 1:n
                for l = 1:16
                    c(q)   = tessellation{i,j}(l,1)^2/r^2 + tessellation{i,j}(l,2)^2/(1/4*r^2) - 1;
                    c(q+1) = abs(tessellation{i,j}(l,1)) - HALF_WIDTH;
                    q = q + 2;
                end
            end
        end
        ceq = zeros(1, 4*m + 4*n);
        k = 1;
        for i = 1:m      % left and right edges clamped straight
            ceq(k)   = tessellation{i,1}(14,1) + HALF_WIDTH;
            ceq(k+1) = tessellation{i,1}(9,1)  + HALF_WIDTH;
            ceq(k+2) = tessellation{i,n}(3,1)  - HALF_WIDTH;
            ceq(k+3) = tessellation{i,n}(8,1)  - HALF_WIDTH;
            k = k + 4;
        end
        for j = 1:n      % bottom and top follow the ellipse
            ceq(k)   = shape(s, tessellation{1,j}(12,:), r);
            ceq(k+1) = shape(s, tessellation{1,j}(5,:),  r);
            ceq(k+2) = shape(s, tessellation{m,j}(15,:), r);
            ceq(k+3) = shape(s, tessellation{m,j}(2,:),  r);
            k = k + 4;
        end

    case 3   % ---- vase ---------------------------------------------------
        ceq = zeros(1, 4*m + 4*n);
        k = 1;
        for i = 1:m      % left and right edges clamped straight
            ceq(k)   = tessellation{i,1}(14,1) + HALF_WIDTH;
            ceq(k+1) = tessellation{i,1}(9,1)  + HALF_WIDTH;
            ceq(k+2) = tessellation{i,n}(3,1)  - HALF_WIDTH;
            ceq(k+3) = tessellation{i,n}(8,1)  - HALF_WIDTH;
            k = k + 4;
        end
        if strcmp(mode, 'nonrigid')
            % explicit form, as written in nonrigid.m (top row first)
            for j = 1:n
                ceq(k)   = -sqrt(1/3*r^2 - 1/3*tessellation{m,j}(15,1)^2) + r - tessellation{m,j}(15,2);
                ceq(k+1) = -sqrt(1/3*r^2 - 1/3*tessellation{m,j}(2,1)^2)  + r - tessellation{m,j}(2,2);
                ceq(k+2) =  sqrt(1/3*r^2 - 1/3*tessellation{1,j}(12,1)^2) - r - tessellation{1,j}(12,2);
                ceq(k+3) =  sqrt(1/3*r^2 - 1/3*tessellation{1,j}(5,1)^2)  - r - tessellation{1,j}(5,2);
                k = k + 4;
            end
        else
            % implicit form via shape.m, as written in rigid.m; the bottom
            % arc is the same curve mirrored, hence the -r
            for j = 1:n
                ceq(k)   = shape(s, tessellation{1,j}(12,:), -r);
                ceq(k+1) = shape(s, tessellation{1,j}(5,:),  -r);
                ceq(k+2) = shape(s, tessellation{m,j}(15,:),  r);
                ceq(k+3) = shape(s, tessellation{m,j}(2,:),   r);
                k = k + 4;
            end
        end

    case 4   % ---- wavy ---------------------------------------------------
        % NOTE the cosine written here is the one that actually governs the
        % result. shape.m's s == 4 branch uses different coefficients and is
        % not used for this shape - see README.
        ceq = zeros(1, 4*n + 4*m);
        k = 1;
        for j = 1:n      % top and bottom follow the cosine
            ceq(k)   = tessellation{1,j}(12,2) + 0.45*cos(0.8*pi*(tessellation{1,j}(12,1)-(pi/(0.8*pi)))) + 1.5;
            ceq(k+1) = tessellation{1,j}(5,2)  + 0.45*cos(0.8*pi*(tessellation{1,j}(5,1)-(pi/(0.8*pi))))  + 1.5;
            ceq(k+2) = tessellation{m,j}(15,2) - 0.45*cos(0.8*pi*(tessellation{m,j}(15,1)-(pi/(0.8*pi)))) - 1.5;
            ceq(k+3) = tessellation{m,j}(2,2)  - 0.45*cos(0.8*pi*(tessellation{m,j}(2,1)-(pi/(0.8*pi))))  - 1.5;
            k = k + 4;
        end
        for i = 1:m      % left and right edges clamped straight
            ceq(k)   = tessellation{i,1}(14,1) + HALF_WIDTH;
            ceq(k+1) = tessellation{i,1}(9,1)  + HALF_WIDTH;
            ceq(k+2) = tessellation{i,n}(3,1)  - HALF_WIDTH;
            ceq(k+3) = tessellation{i,n}(8,1)  - HALF_WIDTH;
            k = k + 4;
        end

    case 5   % ---- heart --------------------------------------------------
        ceq = zeros(1, 4*m + 4*n);
        k = 1;
        for i = 1:m
            ceq(k)   = shape(s, tessellation{i,1}(14,:), r);
            ceq(k+1) = shape(s, tessellation{i,1}(9,:),  r);
            ceq(k+2) = shape(s, tessellation{i,n}(3,:),  r);
            ceq(k+3) = shape(s, tessellation{i,n}(8,:),  r);
            k = k + 4;
        end
        for j = 1:n
            ceq(k)   = shape(s, tessellation{1,j}(12,:), r);
            ceq(k+1) = shape(s, tessellation{1,j}(5,:),  r);
            ceq(k+2) = shape(s, tessellation{m,j}(15,:), r);
            ceq(k+3) = shape(s, tessellation{m,j}(2,:),  r);
            k = k + 4;
        end
end
end

% -------------------------------------------------------------------------
function [c, ceq] = general_target(T, target, m, n)
% Boundary conditions for an arbitrary target shape.
%   equality   : boundary nodes sit on the curve            [Eq. 7]
%   inequality : every node stays inside the curve
% The boundary nodes of each unit are, per create_unit.m:
%   left 14, 9   right 3, 8   bottom 12, 5   top 15, 2

% --- every node inside the target ---------------------------------------
allXY = zeros(m*n*16, 2);
k = 0;
for i = 1:m
    for j = 1:n
        allXY(k+1:k+16, :) = T{i,j};
        k = k + 16;
    end
end
c = target_value(target, allXY).';

% --- boundary nodes on the target ---------------------------------------
clamped = ~isempty(target.clampx);

leftXY  = zeros(2*m,2);  rightXY = zeros(2*m,2);
for i = 1:m
    leftXY(2*i-1,:)  = T{i,1}(14,:);   leftXY(2*i,:)  = T{i,1}(9,:);
    rightXY(2*i-1,:) = T{i,n}(3,:);    rightXY(2*i,:) = T{i,n}(8,:);
end
botXY = zeros(2*n,2);   topXY = zeros(2*n,2);
for j = 1:n
    botXY(2*j-1,:) = T{1,j}(12,:);     botXY(2*j,:) = T{1,j}(5,:);
    topXY(2*j-1,:) = T{m,j}(15,:);     topXY(2*j,:) = T{m,j}(2,:);
end

if clamped
    % left and right are the loading grips: held straight, not on the curve
    w = target.clampx;
    ceq = [ (leftXY(:,1) + w).' , (rightXY(:,1) - w).' , ...
            target_value(target, botXY).' , target_value(target, topXY).' ];
else
    ceq = [ target_value(target, leftXY).' , target_value(target, rightXY).' , ...
            target_value(target, botXY).'  , target_value(target, topXY).' ];
end
end
