function tessellation_transformed = tessellation_deployment(m,n,unit_length,opening_angle)
%TESSELLATION_DEPLOYMENT  Uniformly deploy a rotating-squares kirigami tessellation.
%
%   tessellation_transformed = tessellation_deployment(m,n,unit_length,opening_angle)
%
%   Starts from the compact m-by-n grid built by create_tessellation and
%   rotates the four rigid quadrant squares of every unit about their outer
%   hinge corners, alternating the sense of rotation (-,+,-,+) as in the
%   classical rotating-squares mechanism. The result is centred on the origin.
%
%   INPUTS
%     m, n          : number of units along y and x
%     unit_length   : side length of one unit
%     opening_angle : rotation angle xi of the squares, in radians.
%                     0 gives the compact (as-cut) pattern, pi/4 the fully
%                     open one. This is the single degree of freedom of the
%                     uniform (periodic) deployment.
%
%   OUTPUT
%     tessellation_transformed : m-by-n cell array of 16-by-2 node lists
%                                (node numbering: see create_unit.m)
%
%   Used both to generate the compact pattern for export (opening_angle = 0)
%   and to generate the initial guess handed to tessellation_optimization.

tessellation = create_tessellation(unit_length,m,n); % compact squares tessellation
% Rotate each square in a single unit
[m,n] = size(tessellation);
tessellation_transformed = cell(m,n);
for i = 1:m
    for j = 1:n
        % Opening the squares makes every unit grow by the same amount, so
        % unit (i,j) has to be pushed outwards by a multiple of that growth
        % to keep the grid connected.
        t_grid = [unit_length/2*(sin(opening_angle)+cos(opening_angle)-1)*(j-1)*2, ...
                  unit_length/2*(sin(opening_angle)+cos(opening_angle)-1)*(i-1)*2];
        % translation that keeps the left/right halves of a unit in contact
        t_unit = [unit_length/2*(1-cos(opening_angle)),0];

        square1 = tessellation{i,j}([1,2,3,4],:);      % Q1 top-right
        square2 = tessellation{i,j}([5,6,7,8],:);      % Q2 bottom-right
        square3 = tessellation{i,j}([9,10,11,12],:);   % Q3 bottom-left
        square4 = tessellation{i,j}([13,14,15,16],:);  % Q4 top-left

        % Rotate each quadrant about its own hinge corner, alternating sense:
        % node 4 / 7 (right-mid) and node 10 / 13 (left-mid) are the pivots.
        square_transformed1 = transform_square(square1, -opening_angle, -t_unit+t_grid, tessellation{i,j}(4,:));
        square_transformed2 = transform_square(square2,  opening_angle, -t_unit+t_grid, tessellation{i,j}(7,:));
        square_transformed3 = transform_square(square3, -opening_angle,  t_unit+t_grid, tessellation{i,j}(10,:));
        square_transformed4 = transform_square(square4,  opening_angle,  t_unit+t_grid, tessellation{i,j}(13,:));

        tessellation_transformed{i,j} = [square_transformed1(1:4,:);
            square_transformed2(1:4,:);
            square_transformed3(1:4,:);
            square_transformed4(1:4,:)];
    end
end
% move tessellation to centre (0,0): x from the bottom-left corner of unit
% (1,1) to the bottom-right corner of unit (1,n), y likewise over the rows
xc = (tessellation_transformed{1,1}(9,1) + tessellation_transformed{1,n}(8,1))/2;
yc = (tessellation_transformed{1,1}(12,2) + tessellation_transformed{m,1}(15,2))/2;
for i = 1:m
    for j = 1:n
        tessellation_transformed{i,j} = [tessellation_transformed{i,j}(:,1) - xc, ...
                                         tessellation_transformed{i,j}(:,2) - yc];
    end
end
end
