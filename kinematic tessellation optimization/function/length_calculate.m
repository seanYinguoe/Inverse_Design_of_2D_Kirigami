function edge_length = length_calculate(nodes)
%LENGTH_CALCULATE  Distance between two nodes.
%
%   edge_length = length_calculate(nodes)
%
%   INPUT
%     nodes : 2-by-2, the (x,y) coordinates of two nodes [n1; n2]
%
%   OUTPUT
%     edge_length : Euclidean distance |n2 - n1|
%
%   Building block of the edge conditions Eqs. (4) and (9) of the paper
%   (corresponding edges must have equal length).

edge_length = norm(nodes(2,:) - nodes(1,:));
end
