function edge_diff = edge_diff(nodes,m,n)
%EDGE_DIFF  Objective term: how much the edge lengths vary between units.
%
%   edge_diff = edge_diff(nodes,m,n)
%
%   Edge-length counterpart of angle_diff: sums the squared difference of
%   every quadrant-square edge of unit (i,j) against the same edge of its
%   upper neighbour (i+1,j) and its right neighbour (i,j+1). Zero for a
%   perfectly periodic pattern, so minimising it keeps the panels as close to
%   equally sized as the boundary condition allows.
%
%   INPUTS
%     nodes : (m*n*16)-by-2 flat node array (see nodes_to_units)
%     m, n  : number of units along y and x
%
%   OUTPUT
%     edge_diff : scalar, sum of squared length differences (length^2)
%
%   COVERAGE  As in angle_diff, every adjacent pair is counted exactly once;
%   the previous loop bounds left most of the last row and column out.

tessellation = nodes_to_units(nodes,m,n);

% the 4 edges of each of the 4 quadrant squares, as node pairs
index = [1 2;2 3;3 4;4 1;5 6;6 7;7 8;8 5;9 10;10 11;11 12;12 9;13 14;14 15;15 16;16 13];
length_total = 0;

for i = 1:m
    for j = 1:n
        for edge = 1:16
            l1 = length_calculate(tessellation{i,j}([index(edge,:)],:));
            if i < m    % compare with the unit above
                l2 = length_calculate(tessellation{i+1,j}([index(edge,:)],:));
                length_total = length_total + (l1 - l2)^2;
            end
            if j < n    % compare with the unit to the right
                l3 = length_calculate(tessellation{i,j+1}([index(edge,:)],:));
                length_total = length_total + (l1 - l3)^2;
            end
        end
    end
end
edge_diff = length_total;
end
