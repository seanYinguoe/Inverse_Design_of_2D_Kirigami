function nodes = units_to_nodes(tessellation)
%UNITS_TO_NODES  Flatten the unit cell array into the optimisation variable.
%
%   nodes = units_to_nodes(tessellation)
%
%   Inverse of nodes_to_units: stacks the m-by-n cell array of 16-by-2 units
%   row-major into the single (m*n*16)-by-2 array fmincon optimises over.
%
%   INPUT
%     tessellation : m-by-n cell array of 16-by-2 node lists
%
%   OUTPUT
%     nodes : (m*n*16)-by-2 array; unit (i,j) occupies rows
%             (i-1)*n*16 + (j-1)*16 + (1:16)

[m, n] = size(tessellation);
nodes = zeros(m*n*16,2);
for i = 1:m
    for j = 1:n
        nodes((i-1)*n*16+16*(j-1)+1:(i-1)*n*16+16*j,:) = tessellation{i,j};
    end
end
end
