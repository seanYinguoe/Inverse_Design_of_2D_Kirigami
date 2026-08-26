function tessellation_optimized = tessellation_optimization(tessellation_transformed,s,r,p)
%TESSELLATION_OPTIMIZATION  Stage-1 kinematic optimisation of the cut pattern.
%
%   tessellation_optimized = tessellation_optimization(tessellation_transformed,s,r,p)
%
%   Moves every node of the tessellation with fmincon until the deployed
%   boundary sits on the target shape, subject to the geometric conditions of
%   Section 2.1 of the paper. The design variables are the coordinates of all
%   m*n*16 nodes; the hinges are imposed as LINEAR equalities (Aeq) and the
%   deployability / non-overlap / boundary conditions as NONLINEAR ones
%   (rigid.m or nonrigid.m).
%
%   INPUTS
%     tessellation_transformed : m-by-n cell array, the uniformly deployed
%                                pattern used as the initial guess
%                                (from tessellation_deployment)
%     s : target shape selector, see shape.m (1 circle ... 5 heart)
%     r : characteristic size of the target shape
%     p : deployability mode, 1 = rigid deployable, 2 = non-rigid deployable
%
%   OUTPUT
%     tessellation_optimized : m-by-n cell array of 16-by-2 node lists
%
%   Feed the result to tessellation_compaction to obtain the as-cut pattern.

%% Define the objective function
unit_length = 1;
[m , n] = size(tessellation_transformed);
% number of interior corners the two objective terms are summed over,
% used to normalise the objective so it does not scale with the grid size
M = (m-1)*(n-1)*2 + (m-1)+(n-1);
% keep the pattern as close to periodic as the boundary condition allows:
% angle_diff has units of rad^2 and is converted to a length by unit_length/pi
fun = @(nodes) 1/M * (unit_length/pi * (angle_diff(nodes,m,n)) + edge_diff(nodes,m,n));

%% Define the initial configuration
% fmincon optimises one flat (m*n*16)-by-2 array of node coordinates
tessellation_initial = units_to_nodes(tessellation_transformed);

%% Define the nonlinear conditions (deployability, non-overlap, boundary)
if p == 1
    nonlcon = @(nodes) rigid(nodes,s,m,n,r);
elseif p == 2
    nonlcon = @(nodes) nonrigid(nodes,s,m,n,r);
end
A = []; % linear inequality constraints
b = []; % linear inequality constraints
lb = []; % Lower bounds
ub = []; % Upper bounds

%% Hinge conditions as linear equalities: Aeq*x = 0
% fmincon flattens the (m*n*16)-by-2 variable column-major, so the unknown
% vector is [all x coordinates; all y coordinates] and the y coordinate of
% node k lives at column k + m*n*16. Each hinge contributes two rows
% (x and y): "these two duplicated nodes must stay at the same place".
% Rows beyond the ones filled below stay all-zero and are trivially satisfied.
Aeq = zeros(2*m*n*16,2*m*n*16); % linear equality constraints
beq = zeros(2*m*n*16,1); % linear equality constraints
% hinge node pairs (see the node map in create_unit.m)
%   k = 1..4 : hinges inside a unit   4-7 right-mid, 10-13 left-mid,
%                                     6-11 and 16-1 at the centre
%   k = 5..6 : hinges to the unit on the right  8-9 and 3-14
%   k = 7..8 : hinges to the unit above        15-12 and 2-5
index1 = [4 10 6 16 8 3 15 2];
index2 = [7 13 11 1 9 14 12 5];
eq_row = 1;
% hinges internal to each unit
for i = 1:m
    for j = 1:n
        for k = 1:4
            Aeq(eq_row,index1(k)+(i-1)*n*16+16*(j-1)) = 1;  % x cooridinates
            Aeq(eq_row,index2(k)+(i-1)*n*16+16*(j-1)) = -1;
            Aeq(eq_row+1,index1(k)+(i-1)*n*16+16*(j-1)+m*n*16) = 1; % y coordinates
            Aeq(eq_row+1,index2(k)+(i-1)*n*16+16*(j-1)+m*n*16) = -1;
            eq_row = eq_row+2;
        end
    end
end
% hinges between horizontally adjacent units (j and j+1)
for i = 1:m
    for j = 1:n-1
        for k = 5:6
            Aeq(eq_row,index1(k)+(i-1)*n*16+16*(j-1)) = 1;
            Aeq(eq_row,index2(k)+(i-1)*n*16+16*j) = -1;
            Aeq(eq_row+1,index1(k)+(i-1)*n*16+16*(j-1)+m*n*16) = 1;
            Aeq(eq_row+1,index2(k)+(i-1)*n*16+16*j+m*n*16) = -1;
            eq_row = eq_row+2;
        end
    end
end
% hinges between vertically adjacent units (i and i+1)
for i = 1:m-1
    for j = 1:n
        for k = 7:8
            Aeq(eq_row,index1(k)+(i-1)*n*16+16*(j-1)) = 1;
            Aeq(eq_row,index2(k)+i*n*16+16*(j-1)) = -1;
            Aeq(eq_row+1,index1(k)+(i-1)*n*16+16*(j-1)+m*n*16) = 1;
            Aeq(eq_row+1,index2(k)+i*n*16+16*(j-1)+m*n*16) = -1;
            eq_row = eq_row+2;
        end
    end
end

%% Solve
% Define the options for the solver
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point','MaxFunEvals',8000);
%options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point',  'StepTolerance', 1e-6);
% Call the fmincon function to minimize the objective function
[nodes_optimized, fval] = fmincon(fun, tessellation_initial, A, b, Aeq, beq, lb, ub, nonlcon, options);
fprintf('Final objective value: %.6g\n', fval);

% back from the flat node array to the m-by-n cell array of units
tessellation_optimized = nodes_to_units(nodes_optimized,m,n);
end
