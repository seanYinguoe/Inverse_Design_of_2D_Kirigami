function points = create_unit(vertices)
%CREATE_UNIT  Split one square into the 4 rigid quadrant squares of a kirigami unit.
%
%   points = create_unit(vertices)
%
%   INPUT
%     vertices : 4-by-2, the corners of the unit square in counter-clockwise
%                order starting bottom-left:
%                   v1 = bottom-left, v2 = bottom-right,
%                   v3 = top-right,   v4 = top-left.
%
%   OUTPUT
%     points   : 16-by-2 node list. Rows 1:4, 5:8, 9:12 and 13:16 are the
%                four rigid quadrant squares, each listed clockwise. This
%                16-node layout is THE data structure of the whole program:
%                every index list in the objective, constraint and plotting
%                functions refers to it.
%
%   NODE MAP OF ONE UNIT  (Q1..Q4 are the four rigid quadrant squares)
%
%          14 --------15/2-------- 3          v4 = 14 (top-left)
%           |    Q4    |    Q1     |          v3 =  3 (top-right)
%           |          |           |          v2 =  8 (bottom-right)
%          13        16/1          4          v1 =  9 (bottom-left)
%          /10        11/6         7\         (mid-edge and centre nodes are
%           |    Q3    |    Q2     |           duplicated, one copy per square)
%           |          |           |
%           9 --------12/5-------- 8
%
%     Q1 = rows 1:4   = [centre, top-mid, top-right corner, right-mid]
%     Q2 = rows 5:8   = [bottom-mid, centre, right-mid, bottom-right corner]
%     Q3 = rows 9:12  = [bottom-left corner, left-mid, centre, bottom-mid]
%     Q4 = rows 13:16 = [left-mid, top-left corner, top-mid, centre]
%
%   DUPLICATED NODES AND HINGES
%     A node shared by two squares appears twice, once in each square. Whether
%     the two copies are tied together (a hinge / ligament) or free to separate
%     (a cut) is what makes the pattern a kirigami rather than a solid sheet:
%
%       tied within a unit (hinges):    4<->7  (right-mid),  10<->13 (left-mid)
%                                       1<->16 (centre),      6<->11 (centre)
%       tied across units (hinges):      8<->9,  3<->14  (unit j to unit j+1)
%                                       15<->12, 2<->5   (unit i to unit i+1)
%       never tied (these open up):     2<->15 and 5<->12 inside a unit, and
%                                       the two centre pairs relative to each
%                                       other - this is the rotating-squares
%                                       void that lets the pattern deploy.
%
%     The same tie list is imposed as linear equalities in
%     tessellation_optimization.m (index1 / index2).

% Midpoints of the four sides: mid(k) is the midpoint of edge v(k)->v(k+1),
% i.e. mid1 = bottom, mid2 = right, mid3 = top, mid4 = left.
midpoints = (vertices + [vertices(2:end,:); vertices(1,:)])/2;

% Centre of the unit, shared (as 4 copies) by all quadrant squares.
center = mean(vertices);

% Assemble the 4 quadrant squares, each listed clockwise (see node map above).
points = [center;midpoints(3,:);vertices(3,:);midpoints(2,:);...     % Q1 top-right
    midpoints(1,:);center;midpoints(2,:);vertices(2,:);...           % Q2 bottom-right
    vertices(1,:);midpoints(4,:);center;midpoints(1,:);...           % Q3 bottom-left
    midpoints(4,:);vertices(4,:);midpoints(3,:);center];             % Q4 top-left
end
