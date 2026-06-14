function length = length_calculate(nodes)
% node1, node2are 2x1 matrices representing the (x,y) coordinates of each node

% Calculate the vectors formed by the nodes
vector = nodes(2,:) - nodes(1,:);

% Calculate the magnitudes of the vectors
length = norm(vector);

end
