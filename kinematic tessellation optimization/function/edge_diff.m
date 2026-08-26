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

tessellation = nodes_to_units(nodes,m,n);

% the 4 edges of each of the 4 quadrant squares, as node pairs
index = [1 2;2 3;3 4;4 1;5 6;6 7;7 8;8 5;9 10;10 11;11 12;12 9;13 14;14 15;15 16;16 13];
length_total = 0;

% compare every unit with the one above and the one to its right
for i = 1:m-1
    for j = 1:n-1
        for edge = 1:16
            nodes1 = tessellation{i,j}([index(edge,:)],:);
            nodes2 = tessellation{i+1,j}([index(edge,:)],:);
            nodes3 = tessellation{i,j+1}([index(edge,:)],:);
            length1 = length_calculate(nodes1);
            length2 = length_calculate(nodes2);   % upper neighbour
            length3 = length_calculate(nodes3);   % right neighbour
            length_total = length_total + (length1 - length2)^2 + (length1 - length3)^2;
        end
    end
end
% the loop above skips the last row and column; close the sum at the far
% corner unit (m,n) by comparing it with its lower and left neighbours
for edge = 1:16
    nodes1 = tessellation{m,n}([index(edge,:)],:);
    nodes2 = tessellation{m-1,n}([index(edge,:)],:);
    nodes3 = tessellation{m,n-1}([index(edge,:)],:);
    length1 = length_calculate(nodes1);
    length2 = length_calculate(nodes2);
    length3 = length_calculate(nodes3);
    length_total = length_total + (length1 - length2)^2 + (length1 - length3)^2;
end
edge_diff = length_total;
end
