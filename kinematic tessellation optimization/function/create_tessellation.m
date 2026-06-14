function tessellation = create_tessellation(l,m,n)
%   generates n by n squares tessellation
%   l: the length of each square
%   n: the number of squares in each row/column
%   tessellation: a cell structure that contains the information of each
%                 kirigami units and cooresponding points

% initialize the cell array to store the nodes of each square
tessellation = cell(m, n);

% loop through each row and column to generate the nodes of each square
for i = 1:m
    for j = 1:n
        % calculate the x and y coordinates of the four vertices of the square
        x = (j-1)*l + [0 l l 0];
        y = (i-1)*l + [0 0 l l];
        
        % store the coordinates in a 4-by-2 matrix
        nodes = [x' y'];
        
        % store the matrix in the corresponding cell of the cell array
        tessellation{i,j} = create_unit(nodes);
    end
end
end
