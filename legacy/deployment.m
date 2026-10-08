% DEPLOYMENT  Scratch script: replay the deployment of a single unit.
%
%   SCRIPT (not a function) - run it after main.m has produced
%   tessellation_optimized and tessellation_compacted in the base workspace.
%
%   Takes unit (1,1) in its compacted state and rotates its four quadrant
%   squares back to the optimised deployed position, then plots the result.
%   Used to check the deployment path of one unit by hand; the loop bounds
%   are hard-coded to 1:1, widen them to animate more units.
%
%   Node numbering: see create_unit.m.

tessellation = tessellation_optimized;
tessellation_transformed = tessellation_optimized;
% rotate each quadrant square from the compacted state to the deployed one
for i = 1:1
    for j = 1:1
        d1 = tessellation_compacted{i,j}(4,:) - tessellation{i,j}(4,:);
        d2 = tessellation_compacted{i,j}(13,:) - tessellation{i,j}(13,:);
        angle1 = angle_calculate([tessellation_optimized{i,j}(1,:)+d1;tessellation_optimized{i,j}(4,:)+d1;tessellation_compacted{i,j}(1,:)]);
        angle2 = angle_calculate([tessellation_optimized{i,j}(6,:)+d1;tessellation_optimized{i,j}(7,:)+d1;tessellation_compacted{i,j}(6,:)]);
        angle3 = angle_calculate([tessellation_optimized{i,j}(11,:)+d2;tessellation_optimized{i,j}(10,:)+d2;tessellation_compacted{i,j}(11,:)]);
        angle4 = angle_calculate([tessellation_optimized{i,j}(16,:)+d2;tessellation_optimized{i,j}(13,:)+d2;tessellation_compacted{i,j}(16,:)]);
        square_transformed1 = transform_square(tessellation{i,j}(1:4,:),angle1,d1,tessellation{i,j}(4,:));
        square_transformed2 = transform_square(tessellation{i,j}(5:8,:),-angle2,d1,tessellation{i,j}(7,:));
        square_transformed3 = transform_square(tessellation{i,j}(9:12,:),angle3,d2,tessellation{i,j}(10,:));
        square_transformed4 = transform_square(tessellation{i,j}(13:16,:),-angle4,d2,tessellation{i,j}(13,:));
        tessellation_transformed{i,j} = [square_transformed1(1:4,:);
            square_transformed2(1:4,:);
            square_transformed3(1:4,:);
            square_transformed4(1:4,:)];
    end
end
figure(1)
plot_tessellation(tessellation_transformed)
