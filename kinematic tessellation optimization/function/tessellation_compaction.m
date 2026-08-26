function compact_tessellation = tessellation_compaction(tessellation)
%TESSELLATION_COMPACTION  Fold a deployed tessellation back to its as-cut state.
%
%   compact_tessellation = tessellation_compaction(tessellation)
%
%   Inverse of the deployment: takes the optimised deployed configuration and
%   closes every hinge so the panels touch, which recovers the flat cut
%   pattern that has to be laser-cut. Because the optimised pattern is not
%   periodic, the units are closed one at a time and then re-assembled:
%
%     1. close the four quadrant squares of each unit about their hinges;
%     2. place the first row of units left to right, each rotated so it
%        touches its left neighbour;
%     3. place the first column of units bottom to top likewise;
%     4. fill in the remaining units, each keyed to the unit on its left and
%        the unit below it;
%     5. re-centre the whole pattern.
%
%   Throughout, angle_calculate returns an unsigned angle, so ifoverlapping
%   is used to decide which way the panel has to be rotated.
%
%   INPUT
%     tessellation : m-by-n cell array of 16-by-2 node lists, deployed state
%                    (typically the output of tessellation_optimization)
%
%   OUTPUT
%     compact_tessellation : m-by-n cell array of 16-by-2 node lists,
%                            compact (as-cut) state
%
%   Node numbering: see create_unit.m.

[m, n] = size(tessellation);
compact_tessellation = cell(m,n);

%% 1. Close the four quadrant squares within each unit
for i = 1:m
    for j = 1:n
        % opening angles at the top-mid (2-1-15) and left-mid (11-10-16) hinges
        angle1 = angle_calculate(tessellation{i,j}([2,1,15],:));
        angle2 = angle_calculate(tessellation{i,j}([11,10,16],:));
        % Q1: rotate about its centre node (node 1) until it closes
        square_transformed1 = transform_square(tessellation{i,j}(1:4,:),angle1,[0 0],tessellation{i,j}(1,:));
        % Q2: first slide it so its hinge node 7 meets node 4 of the closed Q1
        t = square_transformed1(4,:) - tessellation{i,j}(7,:);
        square_transformed2 = transform_square(tessellation{i,j}(5:8,:),0,t,tessellation{i,j}(7,:));
        % then rotate it about that hinge until it closes onto Q1
        nodes3 = [square_transformed1([1,4],:);square_transformed2(2,:)];
        angle3 = angle_calculate(nodes3);
        if ifoverlapping(nodes3) >= 0
            angle3 = -angle3;   % turn the other way round
        end
        square_transformed2 = transform_square(tessellation{i,j}(5:8,:),angle3,t,tessellation{i,j}(7,:));
        % Q4 is the reference panel and is left where it is
        square_transformed4 = tessellation{i,j}(13:16,:);
        % Q3: rotate about its left-mid hinge (node 10)
        square_transformed3 = transform_square(tessellation{i,j}(9:12,:),angle2,[0 0],tessellation{i,j}(10,:));
        compact_tessellation{i,j} = [square_transformed1;square_transformed2;square_transformed3;square_transformed4];
    end
end

%% 2. Place the first row of units, left to right
% put unit (1,1) at [-n/2 -m/2] and align its bottom edge (9->12) with +x
t = [-n/2 -m/2] - compact_tessellation{1,1}(9,:);
nodes = [compact_tessellation{1,1}([12,9],:);1e11,0];   % 1e11 ~ a point far along +x
closing_angle = angle_calculate(nodes);
compact_tessellation{1,1} = transform_square(compact_tessellation{1,1}(1:16,:),-closing_angle,t,compact_tessellation{1,1}(9,:));
for j = 2:n
    % slide unit (1,j) so its bottom-left corner (node 9) meets the
    % bottom-right corner (node 8) of its left neighbour
    t = compact_tessellation{1,j-1}(8,:) - compact_tessellation{1,j}(9,:);
    compact_tessellation{1,j} = transform_square(compact_tessellation{1,j}(1:16,:),0,t,compact_tessellation{1,j}(9,:));
    % then rotate it about that shared corner until the two edges coincide
    nodes = [compact_tessellation{1,j-1}(7,:);compact_tessellation{1,j}([9,10],:)];
    closing_angle = angle_calculate(nodes);
    if ifoverlapping(nodes) >= 0
        closing_angle = -closing_angle;
    end
    compact_tessellation{1,j} = transform_square(compact_tessellation{1,j}(1:16,:),closing_angle,[0 0],compact_tessellation{1,j}(9,:));
end

%% 3. Place the first column of units, bottom to top
for i = 2:m
    % node 9 of unit (i,1) meets node 14 (top-left corner) of the unit below
    t = compact_tessellation{i-1,1}(14,:) - compact_tessellation{i,1}(9,:);
    compact_tessellation{i,1} = transform_square(compact_tessellation{i,1}(1:16,:),0,t,compact_tessellation{i,1}(9,:));
    nodes = [compact_tessellation{i-1,1}(15,:);compact_tessellation{i,1}([9,12],:)];
    closing_angle = angle_calculate(nodes);
    if ifoverlapping(nodes) >= 0
        closing_angle = -closing_angle;
    end
    compact_tessellation{i,1} = transform_square(compact_tessellation{i,1}(1:16,:),closing_angle,[0 0],compact_tessellation{i,1}(9,:));
end

%% 4. Fill in the remaining units
for i = 2:m
    for j = 2:n
        % translate against the left neighbour, rotate against the one below
        t = compact_tessellation{i,j-1}(8,:) - compact_tessellation{i,j}(9,:);
        compact_tessellation{i,j} = transform_square(compact_tessellation{i,j}(1:16,:),0,t,compact_tessellation{i,j}(9,:));
        nodes = [compact_tessellation{i-1,j}(15,:);compact_tessellation{i,j}([9,12],:)];
        closing_angle = angle_calculate(nodes);
        if ifoverlapping(nodes) >= 0
            closing_angle = -closing_angle;
        end
        compact_tessellation{i,j} = transform_square(compact_tessellation{i,j}(1:16,:),closing_angle,[0 0],compact_tessellation{i,j}(9,:));
    end
end

%% 5. Translate the whole tessellation to the centre
length_x = compact_tessellation{1,n}(8,1) - compact_tessellation{1,1}(9,1);
length_y = compact_tessellation{m,1}(14,2) - compact_tessellation{1,1}(9,2);
for i = 1:m
    for j = 1:n
        t = [n/2-length_x/2,m/2-length_y/2];
        compact_tessellation{i,j} = transform_square(compact_tessellation{i,j}(1:16,:),0,t,compact_tessellation{i,j}(9,:));
    end
end
end
