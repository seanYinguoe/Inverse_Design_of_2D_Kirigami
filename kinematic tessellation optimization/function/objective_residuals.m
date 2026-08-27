function r = objective_residuals(nodes, m, n, w_angle, w_edge)
%OBJECTIVE_RESIDUALS  The optimisation objective written as a residual vector.
%
%   r = objective_residuals(nodes, m, n, w_angle, w_edge)
%
%   The objective of tessellation_optimization is
%
%       f = 1/M * ( l/pi * angle_diff + edge_diff )
%
%   and both terms are sums of SQUARED differences between neighbouring units.
%   Writing those differences out as a vector r, weighted so that
%
%       f = sum(r.^2)
%
%   buys the gradient cheaply: df/dx = 2 * J' * r, where J is the Jacobian of
%   r. J is sparse - each residual involves two units - so it can be built by
%   sparse_fd's colouring scheme in a few dozen evaluations instead of
%   nvars+1. Without this the objective's own finite differences would become
%   the bottleneck once the constraint Jacobian is made sparse.
%
%   The residual ordering is angles first, then edges; within each, the loop
%   order matches angle_diff.m and edge_diff.m exactly, so
%   sum(r.^2) reproduces the scalar objective to rounding.
%
%   INPUTS
%     nodes   : (m*n*16)-by-2 flat node array, or its column-vector form
%     m, n    : number of units along y and x
%     w_angle : sqrt(l/(M*pi)) - weight applied to each angle difference
%     w_edge  : sqrt(1/M)      - weight applied to each edge difference
%
%   OUTPUT
%     r : column vector of weighted residuals

if size(nodes,2) ~= 2
    nodes = reshape(nodes, [], 2);
end
T = nodes_to_units(nodes, m, n);

idxA = [1 2 3;2 3 4;3 4 1;4 1 2;5 6 7;6 7 8;7 8 5;8 5 6;9 10 11;10 11 12; ...
        11 12 9;12 9 10;13 14 15;14 15 16;15 16 13;16 13 14];
idxE = [1 2;2 3;3 4;4 1;5 6;6 7;7 8;8 5;9 10;10 11;11 12;12 9; ...
        13 14;14 15;15 16;16 13];

% every adjacent pair contributes one residual per corner / per edge
npair = (m-1)*n + m*(n-1);
r = zeros(npair*16*2, 1);
k = 0;

% --- angle residuals -----------------------------------------------------
for i = 1:m
    for j = 1:n
        for q = 1:16
            a1 = angle_calculate(T{i,j}(idxA(q,:),:));
            if i < m
                k = k + 1;
                r(k) = w_angle * (a1 - angle_calculate(T{i+1,j}(idxA(q,:),:)));
            end
            if j < n
                k = k + 1;
                r(k) = w_angle * (a1 - angle_calculate(T{i,j+1}(idxA(q,:),:)));
            end
        end
    end
end

% --- edge residuals ------------------------------------------------------
for i = 1:m
    for j = 1:n
        for q = 1:16
            l1 = length_calculate(T{i,j}(idxE(q,:),:));
            if i < m
                k = k + 1;
                r(k) = w_edge * (l1 - length_calculate(T{i+1,j}(idxE(q,:),:)));
            end
            if j < n
                k = k + 1;
                r(k) = w_edge * (l1 - length_calculate(T{i,j+1}(idxE(q,:),:)));
            end
        end
    end
end

r = r(1:k);
end
