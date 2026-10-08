function plot_tessellation(tessellation)
%PLOT_TESSELLATION  Draw a kirigami tessellation as filled rigid panels.
%
%   plot_tessellation(tessellation)
%
%   Draws the four rigid quadrant squares of every unit as filled patches on
%   the current axes; the cuts and the rotating-squares voids appear as the
%   uncovered background between them.
%
%   INPUT
%     tessellation : m-by-n cell array of 16-by-2 node lists
%                    (node numbering: see create_unit.m)

[m, n] = size(tessellation);
panel_color = [1 .98 .54];   % pale yellow fill for the rigid panels
for i = 1:m
    for j = 1:n
        vertices1 = tessellation{i,j}([1,2,3,4],:);      % Q1 top-right
        vertices2 = tessellation{i,j}([5,6,7,8],:);      % Q2 bottom-right
        vertices3 = tessellation{i,j}([9,10,11,12],:);   % Q3 bottom-left
        vertices4 = tessellation{i,j}([13,14,15,16],:);  % Q4 top-left
        hold on
        patch(vertices1(:,1),vertices1(:,2),panel_color);
        patch(vertices2(:,1),vertices2(:,2),panel_color);
        patch(vertices3(:,1),vertices3(:,2),panel_color);
        patch(vertices4(:,1),vertices4(:,2),panel_color);
        axis equal
    end
end
end
