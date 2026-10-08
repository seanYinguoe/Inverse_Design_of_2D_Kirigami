function angle_diff = angle_diff(nodes,m,n)
%ANGLE_DIFF  Objective term: how much the corner angles vary between units.
%
%   angle_diff = angle_diff(nodes,m,n)
%
%   Sums the squared difference of every corner angle of unit (i,j) against
%   the same corner of its upper neighbour (i+1,j) and its right neighbour
%   (i,j+1). It is zero for a perfectly periodic pattern and grows as the
%   optimiser has to make the units differ from one another, so minimising it
%   keeps the cut pattern as close to uniform as the boundary condition allows
%   (the 1/M * l/pi weighting is applied in tessellation_optimization).
%
%   INPUTS
%     nodes : (m*n*16)-by-2 flat node array (see nodes_to_units)
%     m, n  : number of units along y and x
%
%   OUTPUT
%     angle_diff : scalar, sum of squared angle differences (rad^2)
%
%   COVERAGE  Every adjacent pair of units is counted exactly once. An
%   earlier version looped i = 1:m-1, j = 1:n-1 and then patched in a single
%   extra term for unit (m,n), which left most of the last row and last
%   column out of the objective entirely - those units were free to distort
%   without penalty, and they are exactly the units carrying the boundary
%   condition. The pattern below has no such gap.

tessellation = nodes_to_units(nodes,m,n);
angle_total = 0;

% the 4 corners of each of the 4 quadrant squares, as node triples
% [previous, corner, next] - the angle is measured at the middle entry
index = [1 2 3;2 3 4;3 4 1;4 1 2;5 6 7;6 7 8;7 8 5;8 5 6;9 10 11;10 11 12;
    11 12 9;12 9 10;13 14 15;14 15 16;15 16 13;16 13 14];

for i = 1:m
    for j = 1:n
        for corner = 1:16
            nodes1 = tessellation{i,j}([index(corner,:)],:);
            a1 = angle_calculate(nodes1);
            if i < m    % compare with the unit above
                a2 = angle_calculate(tessellation{i+1,j}([index(corner,:)],:));
                angle_total = angle_total + (a1 - a2)^2;
            end
            if j < n    % compare with the unit to the right
                a3 = angle_calculate(tessellation{i,j+1}([index(corner,:)],:));
                angle_total = angle_total + (a1 - a3)^2;
            end
        end
    end
end
angle_diff = angle_total;
end
