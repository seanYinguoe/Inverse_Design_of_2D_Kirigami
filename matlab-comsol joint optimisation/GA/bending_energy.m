function U = bending_energy(theta)
% theta
% Program for calculate the bending energy of kirigami tessellation analytically
% input: node coordinates nodes = [x1 x2 x3 x4
%                                  y1 y2 y3 y4]
% output: bending energy
% Calculate the opening angle of structure

% Input material and geometric properties
E = 0.05e9;   % the elastic module of material
t = 0.003;    % the thicknes of kirigami sheet
d = 0.0165;     %the gap between cuts 0.02 0.016
w = 0.02;     %the width of cuts 0.02

U = zeros(1,51);

% Calculate the energy of each step
for i = 1:size(theta,2)
    U(i) = (1/24) * E * t * d^3 * (theta(i))^2 / w;
end

% % Get the node coordinate
% node1 = nodes(1);
% node2 = nodes(2);
% node3 = nodes(3);
% node4 = nodes(4);
% 
% % Bending hinge in units
% v1 = node1{:} - node2{:};
% v2 = node3{:}  - node4{:};
% 
% % Calculate the energy of each step
% for i = 1:size(v2,2)
%     theta = real(acos(dot(v1(:,i),v2(:,i))/(norm(v1(:,i))*norm(v2(:,i)))));
%     U(i) = (1/24) * E * t * d^3 * (theta)^2 / w;
% end
