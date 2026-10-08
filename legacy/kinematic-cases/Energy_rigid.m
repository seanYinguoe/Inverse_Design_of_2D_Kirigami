% ENERGY_RIGID  Analytical bending energy of a rigid kirigami tessellation.
%
%   SCRIPT (not a function) - run it after main.m has produced
%   tessellation_compacted and tessellation_optimized in the base workspace.
%
%   Implements the rotating-spring model of Section 2.2 of the paper: all the
%   deformation is assumed to localise in the thin ligaments, each of which
%   behaves as a rotational spring, so the total elastic energy is
%
%       U = sum over hinges of  E*t*d^3*xi^2 / (24*w)          [Eqs. 11, 14]
%
%   where xi is the opening angle of that hinge, i.e. the change in corner
%   angle between the compacted and the deployed state.
%
%   INPUT  (read from the base workspace)
%     tessellation_compacted : compact state,  m-by-n cell of 16-by-2 nodes
%     tessellation_optimized : deployed state, same layout
%   OUTPUT (printed)
%     total bending angle [rad] and total bending energy [J]
%
%   Node numbering: see create_unit.m.

% Read the input
cs = tessellation_compacted; % get the nodes of compacted state
ds = tessellation_optimized; % get the nodes of deployed state

% Calculate the opening angle of structure
[m,n] = size(cs);
angle_total = 0;
U_total = 0;

% Input material and geometric properties
E = 4.33e9;   % the elastic module of material
t = 0.002;    % the thicknes of kirigami sheet
d = 0.03;     %the gap between cuts 0.02
w = 0.02;     %the width of cuts 0.02

% bending hinge in units
index = [2,1,15;1,4,6;1,10,11;5,6,12];
for i = 1:m
    for j = 1:n
        for l = 1:4
            angle_c = cs{i,j}([index(l,:)],:);
            angle_d = ds{i,j}([index(l,:)],:);
            opening_angle = abs(angle_calculate(angle_c) - angle_calculate(angle_d));
            angle_total = angle_total + opening_angle;
            U_total = U_total + bending_energy(E, t, d, w, opening_angle);
        end
    end
end
% bending hinge in adjacent units(horizontal)
index = [4,3,13;7,8,10];
for i = 1:m
    for j = 1:n-1
        for l = 1:2
            nodec1 = cs{i,j}(index(l,1),:);
            nodec2 = cs{i,j}(index(l,2),:);
            nodec3 = cs{i,j+1}(index(l,3),:);
            angle_c = [nodec1;nodec2;nodec3];
            noded1 = ds{i,j}(index(l,1),:);
            noded2 = ds{i,j}(index(l,2),:);
            noded3 = ds{i,j+1}(index(l,3),:);
            angle_d = [noded1;noded2;noded3];
            opening_angle = abs(angle_calculate(angle_c) - angle_calculate(angle_d));
            angle_total = angle_total+opening_angle;
            U_total = U_total + bending_energy(E, t, d, w, opening_angle);
        end
    end
end
% bending hinge in adjacent units(vertical)
index = [9,12,14;8,5,3];
for i = 1:m-1
    for j = 1:n
        for l = 1:2
            nodec1 = cs{i+1,j}(index(l,1),:);
            nodec2 = cs{i+1,j}(index(l,2),:);
            nodec3 = cs{i,j}(index(l,3),:);
            angle_c = [nodec1;nodec2;nodec3];
% NOTE  the deployed triple below is read entirely from unit (i+1,j) while
% the compacted triple mixes units (i+1,j) and (i,j). This looks like a
% copy-paste slip - the third deployed node is probably meant to come from
% ds{i,j} so that it matches nodec3 - but it is left as-is because changing
% it would change the published energy values. Verify before reusing.
            noded1 = ds{i+1,j}(index(l,1),:);
            noded2 = ds{i+1,j}(index(l,2),:);
            noded3 = ds{i+1,j}(index(l,3),:);
            angle_d = [noded1;noded2;noded3];
            opening_angle = abs(angle_calculate(angle_c) - angle_calculate(angle_d));
            angle_total = angle_total+opening_angle;
            U_total = U_total + bending_energy(E, t, d, w, opening_angle);
        end
    end
end
% Display the result
fprintf('The total bending angle is: %f \n', angle_total);
fprintf('The total bending energy is: %f J\n', U_total);

%% Calculate the bending energy
function U = bending_energy(E, t, d, w, theta)
    U = (1/24) * E * t * d^3 * (theta)^2 / w;
end
