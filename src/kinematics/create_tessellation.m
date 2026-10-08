function tessellation = create_tessellation(l,m,n)
%CREATE_TESSELLATION  Build an m-by-n grid of compact kirigami units.
%
%   tessellation = create_tessellation(l,m,n)
%
%   INPUTS
%     l : side length of one unit square
%     m : number of units along y (rows)
%     n : number of units along x (columns)
%
%   OUTPUT
%     tessellation : m-by-n cell array. tessellation{i,j} is the 16-by-2 node
%                    list of unit (i,j) produced by create_unit (rows 1:4,
%                    5:8, 9:12, 13:16 = the four rigid quadrant squares).
%
%   Unit (i,j) occupies the square [(j-1)l, jl] x [(i-1)l, il], so i indexes
%   rows bottom-to-top and j indexes columns left-to-right.

% initialize the cell array to store the nodes of each square
tessellation = cell(m, n);

% loop through each row and column to generate the nodes of each square
for i = 1:m
    for j = 1:n
        % corners of unit (i,j), counter-clockwise from the bottom-left one
        x = (j-1)*l + [0 l l 0];
        y = (i-1)*l + [0 0 l l];

        % store the coordinates in a 4-by-2 matrix
        nodes = [x' y'];

        % split the square into its 4 rigid quadrant squares (16 nodes)
        tessellation{i,j} = create_unit(nodes);
    end
end
end
