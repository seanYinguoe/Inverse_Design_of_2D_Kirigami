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

tessellation = nodes_to_units(nodes,m,n);
angle_total = 0;

% the 4 corners of each of the 4 quadrant squares, as node triples
% [previous, corner, next] - the angle is measured at the middle entry
index = [1 2 3;2 3 4;3 4 1;4 1 2;5 6 7;6 7 8;7 8 5;8 5 6;9 10 11;10 11 12;
    11 12 9;12 9 10;13 14 15;14 15 16;15 16 13;16 13 14];

% compare every unit with the one above and the one to its right
for i = 1:m-1
    for j = 1:n-1
        for corner = 1:16
            nodes1 = tessellation{i,j}([index(corner,:)],:);
            nodes2 = tessellation{i+1,j}([index(corner,:)],:);
            nodes3 = tessellation{i,j+1}([index(corner,:)],:);
            d_angle_up    = (angle_calculate(nodes1) - angle_calculate(nodes2))^2;
            d_angle_right = (angle_calculate(nodes1) - angle_calculate(nodes3))^2;
            angle_total = angle_total + d_angle_up + d_angle_right;
        end
    end
end
% the loop above skips the last row and column; close the sum at the far
% corner unit (m,n) by comparing it with its lower and left neighbours
for corner = 1:16
    nodes1 = tessellation{m,n}([index(corner,:)],:);
    nodes2 = tessellation{m-1,n}([index(corner,:)],:);
    nodes3 = tessellation{m,n-1}([index(corner,:)],:);
    d_angle_down = (angle_calculate(nodes1) - angle_calculate(nodes2))^2;
    d_angle_left = (angle_calculate(nodes1) - angle_calculate(nodes3))^2;
    angle_total = angle_total + d_angle_down + d_angle_left;
end
angle_diff = angle_total;
end
