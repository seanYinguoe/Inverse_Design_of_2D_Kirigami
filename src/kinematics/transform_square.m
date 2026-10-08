function square_transformed = transform_square(square,theta,t,rotate_point)
%TRANSFORM_SQUARE  Rotate a set of nodes about a pivot, then translate them.
%
%   square_transformed = transform_square(square,theta,t,rotate_point)
%
%   Applies  p' = T * R(theta) * (p - rotate_point) + rotate_point  in
%   homogeneous coordinates, i.e. Eqs. (1)-(2) of the paper.
%
%   INPUTS
%     square       : k-by-2 list of (x,y) nodes (any k: a 4-node quadrant
%                    square, or all 16 nodes of a unit)
%     theta        : rotation angle in radians, positive counter-clockwise
%     t            : [tx ty] translation applied after the rotation
%     rotate_point : [x y] pivot the rotation is taken about
%
%   OUTPUT
%     square_transformed : k-by-2 list of transformed nodes

% Rotation matrix R(theta), Eq. (1)
R = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
% Translation matrix T, Eq. (2)
T = [1 0 t(1); 0 1 t(2); 0 0 1];
[k, ~] = size(square);

% Move the pivot to the origin so the rotation is taken about it
square_centred = [square(:,1) - rotate_point(1), square(:,2) - rotate_point(2)];

% Apply the transformation in homogeneous coordinates
square_transformed = T * R * [square_centred'; ones(1,k)];

% Convert back to (x,y) coordinates and undo the pivot shift
square_transformed = square_transformed(1:2,:)';
square_transformed = [square_transformed(:,1) + rotate_point(1), ...
                      square_transformed(:,2) + rotate_point(2)];
end
