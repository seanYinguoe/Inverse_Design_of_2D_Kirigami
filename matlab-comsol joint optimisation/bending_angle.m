function theta = bending_angle(nodes)
% Program for calculate the bending angle of kirigami tessellation analytically
% input: node coordinates nodes = [x1 x2 x3 x4
%                                  y1 y2 y3 y4]
% output: bending angle
% Calculate the opening angle of structure

% Get the node coordinate
node1 = nodes(1);
node2 = nodes(2);
node3 = nodes(3);
node4 = nodes(4);

% Bending hinge in units
v1 = node1{:} - node2{:};
v2 = node3{:}  - node4{:};

% Calculate the energy of each step
for i = 1:size(v1,2)
    theta(i) = real(acos(dot(v1(:,i),v2(:,i))/(norm(v1(:,i))*norm(v2(:,i)))));
end


