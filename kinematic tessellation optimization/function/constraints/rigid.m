function [c,ceq] = rigid(nodes,s,m,n,r)
%RIGID  Nonlinear constraints for a RIGID-deployable kirigami tessellation.
%
%   [c,ceq] = rigid(nodes,s,m,n,r)
%
%   Constraint function handed to fmincon by tessellation_optimization when
%   p == 1. Rigid deployability means the pattern deploys by pure rotation of
%   the panels about their hinges, with no stretching: by the theorem used in
%   the paper this holds if and only if every cut forms a PARALLELOGRAM void,
%   which is what the angle and edge blocks below enforce.
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
%   brackets; nonrigid.m has the same structure with the weaker conditions):
%     1. angle constraints    - sum of angles = 2*pi around interior corners
%                               [Eq. 5] AND = pi for colinearity [Eq. 10]
%     2. edge constraints     - corresponding edges equal, so the voids are
%                               parallelograms [Eq. 9]
%     3. symmetry conditions  - the design is mirrored left/right and top/
%                               bottom, which quarters the number of free
%                               parameters
%     4. square condition     - the compacted state must close up into an
%                               n-by-m rectangle [Eq. 8]
%     5. non-overlap          - <v1 x v2, n_hat> >= 0 at every corner [Eq. 6]
%     6. boundary condition   - boundary nodes must lie on the target shape
%                               [Eq. 7], interior nodes must stay inside it
%
%   Node numbering: see create_unit.m. Every "index" array below is a list of
%   node triples (for angles) or pairs (for edges) in that numbering.

% define the constraints function
c = []; % nonlinear inequality matrix condition
ceq = []; % nonlinear equality matrix condition
tessellation = nodes_to_units(nodes,m,n);
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
        % rigid deployability [Eq. 10]: pairs of angles across a cut must sum
        % to pi, i.e. the two edges are colinear in the compacted state, so
        % the void between them is a parallelogram
        ceq(k+1) = angle_calculate(tessellation{i,j}(index1(1,:),:)) + angle_calculate(tessellation{i,j}(index1(4,:),:)) - pi;
        ceq(k+2) = angle_calculate(tessellation{i,j}(index1(5,:),:)) + angle_calculate(tessellation{i,j}(index1(6,:),:)) - pi;
        %eq(k+3) = angle_calculate(tessellation{i,j}(index1(7,:),:)) + angle_calculate(tessellation{i,j}(index1(8,:),:)) - pi;
        k = k+3;
    end
end
% rigid deployability [Eq. 10] across vertically adjacent units
index1 = [1 2 3;8 5 6];
for i = 1:m-1
    for j = 1:n
        ceq(k) = angle_calculate(tessellation{i,j}(index1(1,:),:)) + angle_calculate(tessellation{i+1,j}(index1(2,:),:)) - pi;
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
% % edge constraints
% % contractibility constraints of edges in each units(rhombus)
% index2 = [1 2;15 16;1 4;6 7;13 16;10 11;5 6;11 12];
% for i = 1:m
%     for j = 1:n
%         ceq(k) = length_calculate(tessellation{i,j}(index2(1,:),:)) - length_calculate(tessellation{i,j}(index2(2,:),:));
%         ceq(k+1) = length_calculate(tessellation{i,j}(index2(3,:),:)) - length_calculate(tessellation{i,j}(index2(4,:),:));
%         ceq(k+2) = length_calculate(tessellation{i,j}(index2(5,:),:)) - length_calculate(tessellation{i,j}(index2(6,:),:));
%         ceq(k+3) = length_calculate(tessellation{i,j}(index2(7,:),:)) - length_calculate(tessellation{i,j}(index2(8,:),:));
%         ceq(k+4) = length_calculate(tessellation{i,j}(index2(4,:),:)) - length_calculate(tessellation{i,j}(index2(5,:),:));
%         k = k+5;
%     end
% end
% % contractibility constraints of edges in adjacent units(rhombus)
% index2 = [3 4;13 14;7 8;9 10];
% for i = 1:m
%     for j = 1:n-1
%         ceq(k) = length_calculate(tessellation{i,j}(index2(1,:),:)) - length_calculate(tessellation{i,j}(index2(3,:),:));
%         ceq(k+1) = length_calculate(tessellation{i,j}(index2(1,:),:)) - length_calculate(tessellation{i,j+1}(index2(2,:),:));
%         ceq(k+2) = length_calculate(tessellation{i,j}(index2(3,:),:)) - length_calculate(tessellation{i,j+1}(index2(4,:),:));
%         k = k+3;
%     end
% end
% index2 = [15 16;11 12;2 3;5 8;14 15;9 12];
% for i = 1:m-1
%     for j= 1:n
%         ceq(k) = length_calculate(tessellation{i,j}(index2(1,:),:)) - length_calculate(tessellation{i+1,j}(index2(2,:),:));
%         ceq(k+1) = length_calculate(tessellation{i,j}(index2(3,:),:)) - length_calculate(tessellation{i+1,j}(index2(4,:),:));
%         ceq(k+2) = length_calculate(tessellation{i,j}(index2(5,:),:)) - length_calculate(tessellation{i+1,j}(index2(6,:),:));
%         ceq(k+3) = length_calculate(tessellation{i,j}(index2(3,:),:)) - length_calculate(tessellation{i,j}(index2(5,:),:));
%         k = k+4;
%     end
% end
%% 2. Edge constraints (parallelogram voids)
% Corresponding edges either side of a cut must be equal in length [Eq. 9];
% together with the colinearity conditions above this makes every void a
% parallelogram, which is the necessary and sufficient condition for rigid
% deployability. index2 rows are node pairs defining one edge.
index2 = [1 2;15 16;1 4;6 7;13 16;10 11;5 6;11 12];
for i = 1:m
    for j = 1:n
        % the sum of the two edges are the same
        ceq(k) = length_calculate(tessellation{i,j}(index2(5,:),:)) - length_calculate(tessellation{i,j}(index2(4,:),:));
        ceq(k+1) = length_calculate(tessellation{i,j}(index2(3,:),:)) - length_calculate(tessellation{i,j}(index2(6,:),:));
        %ceq(k+2) = length_calculate(tessellation{i,j}(index2(1,:),:)) - length_calculate(tessellation{i,j}(index2(2,:),:));
        %ceq(k+3) = length_calculate(tessellation{i,j}(index2(7,:),:)) - length_calculate(tessellation{i,j}(index2(8,:),:));
        k = k+2;
    end
end
% same parallelogram condition for the voids straddling two units
index2 = [3 4;13 14;7 8;9 10;1 4;6 7;13 16;10 11];
for i = 1:m
    for j = 1:n-1
        ceq(k) = length_calculate(tessellation{i,j}(index2(1,:),:)) - length_calculate(tessellation{i,j+1}(index2(4,:),:));
        ceq(k+1) = length_calculate(tessellation{i,j}(index2(3,:),:)) - length_calculate(tessellation{i,j+1}(index2(2,:),:));
        k = k+2;
    end
end
index2 = [16 15;11 12;6 5;2 1;9 12;14 15;5 8;2 3];
for i = 1:m-1
    for j= 1:n
        ceq(k) = length_calculate(tessellation{i,j}(index2(1,:),:)) - length_calculate(tessellation{i+1,j}(index2(3,:),:));
        ceq(k+1) = length_calculate(tessellation{i,j}(index2(4,:),:)) - length_calculate(tessellation{i+1,j}(index2(2,:),:));
        k = k+2;
    end
end
index2 = [2 3;14 15;5 8;9 12];
for i = 1:m-1
    for j= 1:n-1
        ceq(k) = length_calculate(tessellation{i+1,j}(index2(3,:),:)) - length_calculate(tessellation{i,j+1}(index2(2,:),:));
        ceq(k+1) = length_calculate(tessellation{i,j}(index2(1,:),:)) - length_calculate(tessellation{i+1,j+1}(index2(4,:),:));
        k = k+2;
    end
end
index2 = [12 11;5 6];
for j = 1:n
    ceq(k) = length_calculate(tessellation{1,j}(index2(1,:),:)) - length_calculate(tessellation{1,j}(index2(2,:),:));
    k = k+1;
end
index2 = [15 16;1 2];
for j = 1:n
    ceq(k) = length_calculate(tessellation{m,j}(index2(1,:),:)) - length_calculate(tessellation{m,j}(index2(2,:),:));
    k = k+1;
end
index2 = [14 15;9 12];
for i = 1:m-1
    ceq(k) = length_calculate(tessellation{i,1}(index2(1,:),:)) - length_calculate(tessellation{i+1,1}(index2(2,:),:));
    k = k+1;
end
index2 = [2 3;5 8];
for i = 1:m-1
    ceq(k) = length_calculate(tessellation{i,n}(index2(1,:),:)) - length_calculate(tessellation{i+1,n}(index2(2,:),:));
    k = k+1;
end

%% 3. Symmetry conditions
% The target shapes are symmetric about both axes, so the design is forced to
% be symmetric too. This is not a physical requirement: it removes about
% three quarters of the independent parameters and makes fmincon converge
% much faster. Mirroring about x reverses the node order within a unit
% (node l <-> node 17-l), which is why the x coordinates are added (they
% cancel) while the y coordinates are subtracted (they must match).
for i = 1:m/2  % left and right symmetry
    for j = 1:(n/2)
        for l = 1:16
            ceq(k) = tessellation{i,j}(l,1) + tessellation{i,n-j+1}(17-l,1);
            ceq(k+1) = tessellation{i,j}(l,2) - tessellation{i,n-j+1}(17-l,2);
            k = k+2;
        end
    end
end
% mirroring about y maps the quadrant squares onto each other as
% Q1<->Q2 and Q3<->Q4, i.e. the node pairs listed here
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
ceq(k+1) = d_total1 - n;
ceq(k+2) = d_total3 - m;
%ceq(k+1) = d_total1/d_total3 - n/m;
k = k+3;
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
%% 6. Boundary condition
% The nodes on the outer edge of the sheet must land on the target curve
% [Eq. 7], evaluated through shape.m, while every other node must stay inside
% it. The boundary nodes are, per unit:
%   left  edge -> nodes 14 and 9    right edge -> nodes 3 and 8
%   bottom edge -> nodes 12 and 5   top   edge -> nodes 15 and 2
% For shapes 2-4 the left and right edges are clamped straight at x = -+2.5
% (the loading grips) and only the top and bottom follow the target curve.
if s == 1  % circle boundary
    % inside the boundary shape
    for i = 1:m
        for j = 1:n
            for l = 1:16
                c(num) = length_calculate([tessellation{i,j}(l,:);[0,0]]) - r;
                %c(num+1) = abs(tessellation{i,j}(l,1)) - 2.5;
                num = num + 1;
            end
        end
    end
    % left and right boundary nodes
    for i = 1:m
        ceq(k) = shape(s,[tessellation{i,1}(14,:)],r);
        ceq(k+1) = shape(s,[tessellation{i,1}(9,:)],r);
        ceq(k+2) = shape(s,[tessellation{i,n}(3,:)],r);
        ceq(k+3) = shape(s,[tessellation{i,n}(8,:)],r);
        k = k+4;
    end
    % left and right boundary shape
    for j = 1:n
        ceq(k) = shape(s,[tessellation{1,j}(12,:)],r);
        ceq(k+1) = shape(s,[tessellation{1,j}(5,:)],r);
        ceq(k+2) = shape(s,[tessellation{m,j}(15,:)],r);
        ceq(k+3) = shape(s,[tessellation{m,j}(2,:)],r);
        k = k+4;
    end

elseif s == 2  % ellipse boundary
    % inside the target shape
    for i = 1:m
        for j = 1:n
            for l = 1:16
                c(num) = tessellation{i,j}(l,1)^2/r^2 + tessellation{i,j}(l,2)^2/(1/4*r^2) - 1;
                c(num+1) = abs(tessellation{i,j}(l,1)) - 2.5;
                num = num + 2;
            end
        end
    end
    % left and right boundary nodes
    for i = 1:m
        ceq(k) = tessellation{i,1}(14,1) + 2.5;
        ceq(k+1) = tessellation{i,1}(9,1) + 2.5;
        ceq(k+2) = tessellation{i,n}(3,1) - 2.5;
        ceq(k+3) = tessellation{i,n}(8,1) - 2.5;
        k = k+4;
    end
    % boundary condition for top and bottom condition
    for j = 1:n
        ceq(k) = shape(s,[tessellation{1,j}(12,:)],r);
        ceq(k+1) = shape(s,[tessellation{1,j}(5,:)],r);
        ceq(k+2) = shape(s,[tessellation{m,j}(15,:)],r);
        ceq(k+3) = shape(s,[tessellation{m,j}(2,:)],r);
        k = k+4;
    end

elseif s == 3 % vase boundary
    % left and right boundary nodes
    for i = 1:m
        ceq(k) = tessellation{i,1}(14,1) + 2.5;
        ceq(k+1) = tessellation{i,1}(9,1) + 2.5;
        ceq(k+2) = tessellation{i,n}(3,1) - 2.5;
        ceq(k+3) = tessellation{i,n}(8,1) - 2.5;
        k = k+4;
    end
    % boundary condition for top and bottom condition
%     for j = 1:n
%         ceq(k) = -sqrt(1/3*r^2-1/3*tessellation{m,j}(15,1)^2) + r - tessellation{m,j}(15,2);
%         ceq(k+1) = -sqrt(1/3*r^2-1/3*tessellation{m,j}(2,1)^2) + r - tessellation{m,j}(2,2);
%         ceq(k+2) = sqrt(1/3*r^2-1/3*tessellation{1,j}(12,1)^2) - r - tessellation{1,j}(12,2);
%         ceq(k+3) = sqrt(1/3*r^2-1/3*tessellation{1,j}(5,1)^2) - r - tessellation{1,j}(5,2);
%         k = k+4;
%     end
    for j = 1:n
        ceq(k) = shape(s,[tessellation{1,j}(12,:)],-r);
        ceq(k+1) = shape(s,[tessellation{1,j}(5,:)],-r);
        ceq(k+2) = shape(s,[tessellation{m,j}(15,:)],r);
        ceq(k+3) = shape(s,[tessellation{m,j}(2,:)],r);
        k = k+4;
    end

elseif s == 4 % wavy boundary
    % inside the target shape
%     for j = 1:n
%         for l = 1:16
%             c(num) = -(tessellation{1,j}(l,2) + 0.3*cos(pi*tessellation{1,j}(l,1)) + 1.3);
%             c(num+1) = tessellation{m,j}(l,2) - 0.3*cos(pi*tessellation{m,j}(l,1)) - 1.3;
%             num = num + 2;
%         end
%     end
    % top and bottom boundary conditions
    for j = 1:n
        ceq(k) = tessellation{1,j}(12,2) + 0.45*cos(0.8*pi*(tessellation{1,j}(12,1)-(pi/(0.8*pi)))) + 1.5;
        ceq(k+1) = tessellation{1,j}(5,2) + 0.45*cos(0.8*pi*(tessellation{1,j}(5,1)-(pi/(0.8*pi)))) + 1.5;
        ceq(k+2) = tessellation{m,j}(15,2) - 0.45*cos(0.8*pi*(tessellation{m,j}(15,1)-(pi/(0.8*pi)))) - 1.5;
        ceq(k+3) = tessellation{m,j}(2,2) - 0.45*cos(0.8*pi*(tessellation{m,j}(2,1)-(pi/(0.8*pi)))) - 1.5;
        k = k + 4;
    end
    % left and right boundary shape for wave
    for i = 1:m
        ceq(k) = tessellation{i,1}(14,1) + 2.5;
        ceq(k+1) = tessellation{i,1}(9,1) + 2.5;
        ceq(k+2) = tessellation{i,n}(3,1) - 2.5;
        ceq(k+3) = tessellation{i,n}(8,1) - 2.5;
        k = k+4;
    end
elseif s == 5  % heart shape
    % left and right boundary nodes
    for i = 1:m
        ceq(k) = shape(s,[tessellation{i,1}(14,:)],r);
        ceq(k+1) = shape(s,[tessellation{i,1}(9,:)],r);
        ceq(k+2) = shape(s,[tessellation{i,n}(3,:)],r);
        ceq(k+3) = shape(s,[tessellation{i,n}(8,:)],r);
        k = k+4;
    end
    % left and right boundary shape
    for j = 1:n
        ceq(k) = shape(s,[tessellation{1,j}(12,:)],r);
        ceq(k+1) = shape(s,[tessellation{1,j}(5,:)],r);
        ceq(k+2) = shape(s,[tessellation{m,j}(15,:)],r);
        ceq(k+3) = shape(s,[tessellation{m,j}(2,:)],r);
        k = k+4;
    end
end
end



