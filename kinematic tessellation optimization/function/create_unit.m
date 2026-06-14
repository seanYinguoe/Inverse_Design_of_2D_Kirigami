function points = create_unit(vertices)
% Create a unit square tessellation unit from a square.

% Compute the midpoints of the sides
midpoints = (vertices + [vertices(2:end,:); vertices(1,:)])/2;

% Compute the center point
center = mean(vertices);
% create joints
points = [center;midpoints(3,:);vertices(3,:);midpoints(2,:);...
    midpoints(1,:);center;midpoints(2,:);vertices(2,:);...
    vertices(1,:);midpoints(4,:);center;midpoints(1,:);...
    midpoints(4,:);vertices(4,:);midpoints(3,:);center];
end