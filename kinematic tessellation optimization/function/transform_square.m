function square_transformed = transform_square(square,theta,t,rotate_point)
% input parameters
% square: geometric parameters of square
% theta:Define the rotation angle (in radians)
% l:the length of square
% t=[tx ty]:Define the translation vector
% Create the rotation matrix
R = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
% Create the translation matrix
T = [1 0 t(1); 0 1 t(2); 0 0 1];
[m n] = size(square);
% Translate the point to the origin
suqare_new = [square(:,1) - rotate_point(1),square(:,2) - rotate_point(2)];

% Apply the transformation to the square vertices
square_transformed = T * R * [suqare_new'; ones(1,m)];

% Convert back to (x,y) coordinates
square_transformed = square_transformed(1:2,:)';
square_transformed = [square_transformed(:,1) + rotate_point(1),square_transformed(:,2) + rotate_point(2)];

% plot([square_transformed(:,1);square_transformed(1,1)],[square_transformed(:,2);square_transformed(1,2)])
% Plot the original and transformed squares
% plot(square(:,1), square(:,2), 'b-', square_transformed(:,1), square_transformed(:,2), 'r-');
% axis equal;
end