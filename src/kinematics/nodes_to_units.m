function tessellation = nodes_to_units(nodes,m,n)
%NODES_TO_UNITS  Reshape the flat optimisation vector back into unit cells.
%
%   tessellation = nodes_to_units(nodes,m,n)
%
%   fmincon works on one flat (m*n*16)-by-2 array of coordinates; every
%   geometric function instead wants the m-by-n cell array of 16-by-2 units.
%   This is the inverse of units_to_nodes.
%
%   INPUTS
%     nodes : (m*n*16)-by-2 array, units stacked row-major - unit (i,j)
%             occupies rows (i-1)*n*16 + (j-1)*16 + (1:16)
%     m, n  : number of units along y and x
%
%   OUTPUT
%     tessellation : m-by-n cell array of 16-by-2 node lists
%                    (node numbering: see create_unit.m)

tessellation = cell(m,n);
for i = 1:m
    for j = 1:n
        tessellation{i,j} = nodes((i-1)*n*16+16*(j-1)+1:(i-1)*n*16+16*j,:);
    end
end
end
