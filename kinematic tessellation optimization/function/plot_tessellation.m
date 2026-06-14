function plot_tessellation(tessellation)
% Plots a square tessellation given the tessellation points
[m n] = size(tessellation);
for i = 1:m
    for j = 1:n
        vertices1 = tessellation{i,j}([1,2,3,4],:);
        vertices2 = tessellation{i,j}([5,6,7,8],:);
        vertices3 = tessellation{i,j}([9,10,11,12],:);
        vertices4 = tessellation{i,j}([13,14,15,16],:);
        hold on
        color = [1 .98 .54];
        patch(vertices1(:,1),vertices1(:,2),color);
        patch(vertices2(:,1),vertices2(:,2),color);
        patch(vertices3(:,1),vertices3(:,2),color);
        patch(vertices4(:,1),vertices4(:,2),color);
        axis equal
    end
end
end
        




