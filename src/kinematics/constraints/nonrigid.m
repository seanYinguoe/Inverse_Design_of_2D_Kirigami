function [c,ceq] = nonrigid(nodes,s,m,n,r,opts)
%NONRIGID  Nonlinear constraints for a NON-RIGID deployable kirigami tessellation.
%
%   [c,ceq] = nonrigid(nodes,s,m,n,r)
%
%   Constraint function handed to fmincon by tessellation_optimization when
%   p == 2. Non-rigid deployability only requires that the compacted state
%   closes seamlessly (one intersection per unit); the panels are allowed to
%   be geometrically frustrated on the way, so the voids need NOT be
%   parallelograms. Compared with rigid.m the colinearity conditions
%   (sum of angles = pi, Eq. 10) are dropped and the edge conditions are the
%   weaker Eq. (4) instead of Eq. (9).
%
%   INPUTS
%     nodes : (m*n*16)-by-2 flat node array (see nodes_to_units)
%     s     : target shape selector, see shape.m
%     m, n  : number of units along y and x
%     r     : characteristic size of the target shape
%
%   OUTPUTS
%     c   : inequality constraints, fmincon enforces c   <= 0
%     ceq : equality constraints,   fmincon enforces ceq == 0
%
%   Blocks, in the order they are assembled (paper equation numbers in
%   brackets):
%     1. angle constraints    - sum of angles around every interior corner
%                               equals 2*pi [Eq. 5]
%     2. edge constraints     - corresponding edges have equal length [Eq. 4]
%     3. symmetry conditions  - left/right and top/bottom mirror symmetry
%     4. square condition     - the compacted state closes into an n-by-m
%                               rectangle [Eq. 8]
%     5. non-overlap          - <v1 x v2, n_hat> >= 0 at every corner [Eq. 6]
%     6. boundary condition   - boundary nodes on the target shape [Eq. 7]
%
%
%   OPTIONS (optional 6th argument, struct; omit for the paper's behaviour)
%     MinEdgeLength : lower bound on every quadrant-square edge. The
%                     non-overlap conditions only ask that corners keep their
%                     orientation, which still allows a panel to shrink to
%                     nothing - and it does: on a 4x4 circle the largest and
%                     smallest panels ended up differing by a factor of 400.
%                     A positive value here adds one inequality per edge and
%                     keeps the panels a sensible size. 0 disables it.
%     Symmetry      : when true, add the left-right and top-bottom mirror
%                     conditions. Default FALSE - see the long note at block 3;
%                     switching them on over-constrains the problem and can
%                     make the target boundary unreachable.
%     FreeScale     : when true, the compacted sheet is no longer pinned to
%                     exactly n-by-m. The two absolute-size equalities are
%                     replaced by a single aspect-ratio equality, so the
%                     overall scale becomes a free parameter the optimiser
%                     can use to fit the target boundary. Use together with
%                     fit_initial_guess, which chooses the starting size.
%   Node numbering: see create_unit.m.

% define the constraints function
c = []; % nonlinear inequality matrix condition
ceq = []; % nonlinear equality matrix condition
tessellation = nodes_to_units(nodes,m,n);

% ---- options -------------------------------------------------------------
if nargin < 6 || isempty(opts); opts = struct(); end
if ~isfield(opts,'MinEdgeLength'); opts.MinEdgeLength = 0;     end
if ~isfield(opts,'FreeScale');     opts.FreeScale     = false; end
if ~isfield(opts,'Symmetry');      opts.Symmetry      = false; end

%% 1. Angle constraints
% Contractibility of the angles inside one unit: the four angles meeting at
% the centre of the unit must add up to 2*pi so that the unit closes
% seamlessly when compacted [Eq. 5].
% index1 rows are node triples [previous, corner, next]; the angle is taken
% at the middle node by angle_calculate:
%   rows 1..4 -> the centre corner of Q1, Q2, Q3, Q4
%   rows 5..8 -> the mid-edge corners of Q1, Q2, Q4, Q3
index1 = [4 1 2;5 6 7;10 11 12;15 16 13;3 4 1;6 7 8;16 13 14;9 10 11];
k=1;   % running index into ceq, incremented by every block below
for i = 1:m
    for j = 1:n
        ceq(k) = angle_calculate(tessellation{i,j}(index1(1,:),:)) + angle_calculate(tessellation{i,j}(index1(2,:),:)) + ...
            angle_calculate(tessellation{i,j}(index1(3,:),:)) + angle_calculate(tessellation{i,j}(index1(4,:),:)) - 2*pi;
        k = k+1;
    end
end
% Contractibility of the angles at corners SHARED BY NEIGHBOURING units: the
% angles that meet there must also add up to 2*pi [Eq. 5].
%   index1 -> corner shared by unit (i,j) and its right neighbour  (i,j+1)
%   index2 -> corner shared by unit (i,j) and the unit above it    (i+1,j)
%   index3 -> corner where the four units (i,j), (i+1,j), (i+1,j+1), (i,j+1)
%             meet, one angle contributed by each
index1 = [3 4 1;6 7 8;14 13 16;9 10 11];
index2 = [1 2 3;14 15 16;6 5 8;9 12 11];
index3 = [2 3 4;7 8 5;10 9 12;13 14 15];
for i = 1:m
    for j =1:n-1
        ceq(k) = angle_calculate(tessellation{i,j}(index1(1,:),:)) + angle_calculate(tessellation{i,j}(index1(2,:),:)) + ...
            angle_calculate(tessellation{i,j+1}(index1(3,:),:)) + angle_calculate(tessellation{i,j+1}(index1(4,:),:)) - 2*pi;
        k = k+1;
    end
end
for j = 1:n
    for i = 1:m-1
        ceq(k) = angle_calculate(tessellation{i,j}(index2(1,:),:)) + angle_calculate(tessellation{i,j}(index2(2,:),:)) + ...
            angle_calculate(tessellation{i+1,j}(index2(3,:),:)) + angle_calculate(tessellation{i+1,j}(index2(4,:),:)) - 2*pi;
        k = k+1;
    end
end
for i = 1:m-1
    for j = 1:n-1
        ceq(k) = angle_calculate(tessellation{i,j}(index3(1,:),:)) + angle_calculate(tessellation{i+1,j}(index3(2,:),:)) + ...
            angle_calculate(tessellation{i+1,j+1}(index3(3,:),:)) + angle_calculate(tessellation{i,j+1}(index3(4,:),:)) - 2*pi;
        k = k+1;
    end
end
%% 2. Edge constraints
% Corresponding edges either side of a cut must be equal in length so the cut
% closes seamlessly [Eq. 4]. Unlike rigid.m no colinearity is imposed, so the
% voids may be general quadrilaterals rather than parallelograms.
% index2 rows are node pairs defining one edge.
index2 = [1 2;15 16;1 4;6 7;13 16;10 11;5 6;11 12];
for i = 1:m
    for j = 1:n
        ceq(k) = length_calculate(tessellation{i,j}(index2(3,:),:)) - length_calculate(tessellation{i,j}(index2(4,:),:));
        ceq(k+1) = length_calculate(tessellation{i,j}(index2(1,:),:)) - length_calculate(tessellation{i,j}(index2(2,:),:));
        ceq(k+2) = length_calculate(tessellation{i,j}(index2(5,:),:)) - length_calculate(tessellation{i,j}(index2(6,:),:));
        ceq(k+3) = length_calculate(tessellation{i,j}(index2(7,:),:)) - length_calculate(tessellation{i,j}(index2(8,:),:));
        k = k+4;
    end
end
% same edge condition for the cuts straddling two neighbouring units
index2 = [3 4;13 14;7 8;9 10;1 4;6 7;13 16;10 11];
for i = 1:m
    for j = 1:n-1
        ceq(k) = length_calculate(tessellation{i,j}(index2(1,:),:)) - length_calculate(tessellation{i,j+1}(index2(2,:),:));
        ceq(k+1) = length_calculate(tessellation{i,j}(index2(3,:),:)) - length_calculate(tessellation{i,j+1}(index2(4,:),:));
        k = k+2;
    end
end
index2 = [16 15;11 12;6 5;2 1;9 12;14 15;5 8;2 3];
for i = 1:m-1
    for j= 1:n
        ceq(k) = length_calculate(tessellation{i+1,j}(index2(5,:),:)) - length_calculate(tessellation{i,j}(index2(6,:),:));
        ceq(k+1) = length_calculate(tessellation{i+1,j}(index2(7,:),:)) - length_calculate(tessellation{i,j}(index2(8,:),:));
        k = k+2;
    end
end

%% 3. Symmetry conditions (OPTIONAL - off by default, and for good reason)
% These force the design to be its own mirror image about both axes. They are
% NOT a physical requirement; the paper introduces them only to cut the number
% of independent parameters and speed the solve up.
%
% In practice they over-constrain the problem. Counting the rank of the
% equality Jacobian for a 4x4 circle (historical diagnostic):
%
%            block                     rows   rank added   freedom left
%     ...    edge conditions            ...       ...          115
%            symmetry left-right        128        56           27
%            symmetry top-bottom         64        19            8
%            square conditions           32         8            0
%            BOUNDARY TARGET CURVE       32         0            0   <-- !!
%
% More than half the symmetry rows are redundant with each other, but the
% independent ones still eat the design freedom the boundary condition needs.
% By the time the target-shape equations are added there is nothing left to
% move, so they contribute no rank at all and simply cannot be satisfied.
% That is why the rigid solve used to stall with max|ceq| ~ 1e-1 and produce
% collapsed panels: the system was infeasible, not merely slow.
%
% If you do want a symmetric design, the sound way is to optimise a quarter of
% the sheet and mirror it - reducing the VARIABLES rather than adding
% equations - which is what the paper does for the FEA stage.
if opts.Symmetry
    for i = 1:m/2  % left and right symmetry
        for j = 1:(n/2)
            for l = 1:16
                ceq(k) = tessellation{i,j}(l,1) + tessellation{i,n-j+1}(17-l,1);
                ceq(k+1) = tessellation{i,j}(l,2) - tessellation{i,n-j+1}(17-l,2);
                k = k+2;
            end
        end
    end
    index = [1 6;2 5;3 8;4 7;9 14;10 13;11 16;12 15];
    for i = 1:m/2  % top and bottom symmetry
        for j = 1:(n/2)
            for l = 1:size(index,1)
                ceq(k) = tessellation{i,j}(index(l,1),1) - tessellation{m-i+1,j}(index(l,2),1);
                ceq(k+1) = tessellation{i,j}(index(l,1),2) + tessellation{m-i+1,j}(index(l,2),2);
                k = k+2;
            end
        end
    end
end

%% 4. Square (compact shape) condition
% The compacted state has to close up into an n-by-m rectangle [Eq. 8]:
%   - the summed edge length along the bottom row (d_total1) must equal that
%     along the top row (d_total2), and likewise left (d_total3) vs right
%     (d_total4), otherwise the compacted sheet would not be a rectangle;
%   - d_total1 = n and d_total3 = m fix its actual size;
%   - the "angle of edge" block below forces consecutive boundary edges to be
%     colinear (angle = pi) so each side of the rectangle is straight, and the
%     four corner angles to be pi/2.
% edge
d_total1 = 0;
for j = 1:n
    d = length_calculate(tessellation{1,j}([9,12],:)) + length_calculate(tessellation{1,j}([5,8],:));
    d_total1 = d_total1 + d;
end
d_total2 = 0;
for j = 1:n
    d = length_calculate(tessellation{m,j}([14,15],:)) + length_calculate(tessellation{m,j}([2,3],:));
    d_total2 = d_total2 + d;
end
ceq(k) = d_total1 - d_total2;
k = k+1;
d_total3 = 0;
for i = 1:m
    d = length_calculate(tessellation{i,1}([9,10],:)) + length_calculate(tessellation{i,1}([13,14],:));
    d_total3 = d_total3 + d;
end
d_total4 = 0;
for i = 1:m
    d = length_calculate(tessellation{i,n}([7,8],:)) + length_calculate(tessellation{i,n}([3,4],:));
    d_total4 = d_total4 + d;
end
ceq(k) = d_total3 - d_total4;
if opts.FreeScale
    % Overall size is a free parameter: pin only the aspect ratio, so the
    % optimiser may scale the whole sheet to reach the target boundary.
    % Dropping the two absolute-size rows entirely (rather than setting them
    % to a constant 0) matters - a constant equality contributes an all-zero
    % row to the constraint Jacobian and degrades the KKT conditioning.
    ceq(k+1) = d_total1/n - d_total3/m;
    k = k+2;
else
    % Compacted sheet pinned to exactly n-by-m, as in the paper.
    ceq(k+1) = d_total1 - n;
    ceq(k+2) = d_total3 - m;
    k = k+3;
end
% angle of edge
for j = 1:n
    a = angle_calculate(tessellation{1,j}([9,12,11],:)) + angle_calculate(tessellation{1,j}([6,5,8],:));
    ceq(k) = a - pi;
    k = k+1;
end
for j = 1:n
    a = angle_calculate(tessellation{m,j}([14,15,16],:)) + angle_calculate(tessellation{m,j}([1,2,3],:));
    ceq(k) = a - pi;
    k = k+1;
end
for i = 1:m
    a = angle_calculate(tessellation{i,1}([9,10,11],:)) + angle_calculate(tessellation{i,1}([14,13,16],:));
    ceq(k) = a - pi;
    k = k+1;
end
for i = 1:m
    a = angle_calculate(tessellation{i,n}([3,4,1],:)) + angle_calculate(tessellation{i,n}([8,7,6],:));
    ceq(k) = a - pi;
    k = k+1;
end
for i = 1:m-1
    a = angle_calculate(tessellation{i,1}([13,14,15],:)) + angle_calculate(tessellation{i+1,1}([10,9,12],:));
    ceq(k) = a - pi;
    k = k+1;
end
for i = 1:m-1
    a = angle_calculate(tessellation{i,n}([2,3,4],:)) + angle_calculate(tessellation{i+1,n}([5,8,7],:));
    ceq(k) = a - pi;
    k = k+1;
end
for j = 1:n-1
    a = angle_calculate(tessellation{1,j}([5,8,7],:)) + angle_calculate(tessellation{1,j+1}([12,9,10],:));
    ceq(k) = a - pi;
    k = k+1;
end
for j = 1:n-1
    a = angle_calculate(tessellation{m,j}([2,3,4],:)) + angle_calculate(tessellation{m,j+1}([13,14,15],:));
    ceq(k) = a - pi;
    k = k+1;
end
ceq(k) = angle_calculate(tessellation{1,1}([10,9,12],:)) - pi/2;
ceq(k+1) = angle_calculate(tessellation{m,1}([13,14,15],:)) - pi/2;
ceq(k+2) = angle_calculate(tessellation{1,n}([7,8,5],:)) - pi/2;
ceq(k+3) = angle_calculate(tessellation{m,n}([2,3,4],:)) - pi/2;
k = k+4;
%% 5. Non-overlap constraints
% <v1 x v2, n_hat> >= 0 at every corner between two adjacent panels [Eq. 6],
% written for fmincon as c <= 0 by negating ifoverlapping where the corner is
% expected to turn counter-clockwise. Three groups follow: corners inside a
% unit, corners between horizontally / vertically adjacent units, and finally
% every corner of every quadrant square, which keeps the panels convex
% (a panel that turned inside out would flip the sign).
index3 = [12 6 5;8 7 3;2 16 15;14 13 9;7 6 10;1 7 6;10 1 7;11 10 1];
num = 1;
for i = 1:m
    for j = 1:n
        for l = 1:8
            c(num) = -ifoverlapping(tessellation{i,j}(index3(l,:),:));
            num = num + 1;
        end
    end
end
for i = 1:m
    for j = 1:n-1
        a = tessellation{i,j}(5,:);
        b = tessellation{i,j+1}(12,:);
        f = tessellation{i,j+1}(9,:);
        c(num) = -ifoverlapping([a;f;b]);
        a = tessellation{i,j}(7,:);
        b = tessellation{i,j+1}(10,:);
        f = tessellation{i,j+1}(9,:);
        c(num+1) = ifoverlapping([a;f;b]);
        a = tessellation{i,j}(4,:);
        b = tessellation{i,j+1}(13,:);
        f = tessellation{i,j+1}(14,:);
        c(num+2) = -ifoverlapping([a;f;b]);
        a = tessellation{i,j}(2,:);
        b = tessellation{i,j+1}(15,:);
        f = tessellation{i,j+1}(14,:);
        c(num+3) = ifoverlapping([a;f;b]);
        num = num + 4;
    end
end
for i = 1:m-1
    for j = 1:n
        a = tessellation{i,j}(1,:);
        b = tessellation{i+1,j}(6,:);
        f = tessellation{i+1,j}(5,:);
        c(num) = ifoverlapping([a;f;b]);
        a = tessellation{i,j}(1,:);
        b = tessellation{i+1,j}(11,:);
        f = tessellation{i+1,j}(12,:);
        c(num+1) = -ifoverlapping([a;f;b]);
        a = tessellation{i,j}(3,:);
        b = tessellation{i+1,j}(8,:);
        f = tessellation{i+1,j}(5,:);
        c(num+2) = -ifoverlapping([a;f;b]);
        a = tessellation{i,j}(14,:);
        b = tessellation{i+1,j}(9,:);
        f = tessellation{i+1,j}(12,:);
        c(num+3) = ifoverlapping([a;f;b]);
        num = num + 4;
    end
end
% keep every quadrant square itself convex (all 16 corners turn the same way)
index4 = [1 2 3;2 3 4;3 4 1;4 1 2;5 6 7;6 7 8;7 8 5;8 5 6;9 10 11;10 11 12;11 12 9;12 9 10;13 14 15;14 15 16;15 16 13;16 13 14];
for i = 1:m
    for j = 1:n
        for l = 1:16
            c(num) = -ifoverlapping(tessellation{i,j}(index4(l,:),:));
            num = num + 1;
        end
    end
end
%% mechanic condition
% c(num) = angle_calculate(tessellation{1,1}([16,10,11],:)) -  angle_calculate(tessellation{1,2}([16,10,11],:));
% num = num+1;

%% 6. Boundary condition
% Delegated to boundary_residual so that the constraint files and
% fit_initial_guess.m cannot drift apart. See that function for the list of
% boundary nodes and the per-shape treatment.
[cb, ceqb] = boundary_residual(tessellation, s, r, m, n, 'nonrigid');
if ~isempty(cb)
    c(num:num+numel(cb)-1) = cb;
    num = num + numel(cb);
end
if ~isempty(ceqb)
    ceq(k:k+numel(ceqb)-1) = ceqb;
    k = k + numel(ceqb);
end

%% 7. Minimum edge length (optional)
% The non-overlap conditions above only fix the SENSE of each corner, which
% still permits a panel to shrink towards zero area - and it does in practice.
% This block puts a floor under every quadrant-square edge.
if opts.MinEdgeLength > 0
    index5 = [1 2;2 3;3 4;4 1;5 6;6 7;7 8;8 5;9 10;10 11;11 12;12 9; ...
              13 14;14 15;15 16;16 13];
    for i = 1:m
        for j = 1:n
            for l = 1:16
                c(num) = opts.MinEdgeLength - ...
                    length_calculate(tessellation{i,j}(index5(l,:),:));
                num = num + 1;
            end
        end
    end
end
end
