function compact_tessellation = tessellation_compaction(tessellation)
% get the generalized kirigami pattern from the deployed configuration.
[m n] = size(tessellation);
compact_tessellation = {};
% compact each sqaure units
for i = 1:m
    for j = 1:n
        angle1 = angle_calculate(tessellation{i,j}([2,1,15],:));
        angle2 = angle_calculate(tessellation{i,j}([11,10,16],:));
        square_transformed1 = transform_square(tessellation{i,j}(1:4,:),angle1,[0 0],tessellation{i,j}(1,:));
        t = square_transformed1(4,:) - tessellation{i,j}(7,:);
        square_transformed2 = transform_square(tessellation{i,j}(5:8,:),0,t,tessellation{i,j}(7,:));
        nodes3 = [square_transformed1([1,4],:);square_transformed2(2,:)];
        angle3 = angle_calculate(nodes3);
        if ifoverlapping(nodes3) < 0
            angle3 = angle3;
        else angle3 = -angle3;
        end
        square_transformed2 = transform_square(tessellation{i,j}(5:8,:),angle3,t,tessellation{i,j}(7,:));
        square_transformed4 = tessellation{i,j}(13:16,:);
        square_transformed3 = transform_square(tessellation{i,j}(9:12,:),angle2,[0 0],tessellation{i,j}(10,:));
        compact_tessellation{i,j} = [square_transformed1;square_transformed2;square_transformed3;square_transformed4];
    end
end
% move adjacent square units together(first row)
t = [-n/2 -m/2] - compact_tessellation{1,1}(9,:);
nodes = [compact_tessellation{1,1}([12,9],:);1e11,0];
angle = angle_calculate(nodes);
compact_tessellation{1,1} = transform_square(compact_tessellation{1,1}(1:16,:),-angle,t,compact_tessellation{1,1}(9,:));
for j = 2:n
    t = compact_tessellation{1,j-1}(8,:) - compact_tessellation{1,j}(9,:);
    compact_tessellation{1,j} = transform_square(compact_tessellation{1,j}(1:16,:),0,t,compact_tessellation{1,j}(9,:));
    nodes = [compact_tessellation{1,j-1}(7,:);compact_tessellation{1,j}([9,10],:)];
    angle = angle_calculate(nodes);
    if ifoverlapping(nodes) < 0
        angle = angle;
    else angle = -angle;
    end
    compact_tessellation{1,j} = transform_square(compact_tessellation{1,j}(1:16,:),angle,[0 0],compact_tessellation{1,j}(9,:));
end
% rotate the units(the first column)
for i = 2:m
    t = compact_tessellation{i-1,1}(14,:) - compact_tessellation{i,1}(9,:);
    compact_tessellation{i,1} = transform_square(compact_tessellation{i,1}(1:16,:),0,t,compact_tessellation{i,1}(9,:));
    nodes = [compact_tessellation{i-1,1}(15,:);compact_tessellation{i,1}([9,12],:)];
    angle = angle_calculate(nodes);
    if ifoverlapping(nodes) < 0
        angle = angle;
    else angle = -angle;
    end
    compact_tessellation{i,1} = transform_square(compact_tessellation{i,1}(1:16,:),angle,[0 0],compact_tessellation{i,1}(9,:));
end
% rotate the units(the rest of the units)
for i = 2:m
    for j = 2:n
        t = compact_tessellation{i,j-1}(8,:) - compact_tessellation{i,j}(9,:);
        compact_tessellation{i,j} = transform_square(compact_tessellation{i,j}(1:16,:),0,t,compact_tessellation{i,j}(9,:));
        nodes = [compact_tessellation{i-1,j}(15,:);compact_tessellation{i,j}([9,12],:)];
        angle = angle_calculate(nodes);
        if ifoverlapping(nodes) < 0
            angle = angle;
        else angle = -angle;
        end
        compact_tessellation{i,j} = transform_square(compact_tessellation{i,j}(1:16,:),angle,[0 0],compact_tessellation{i,j}(9,:));
    end
end
% translate the whole tessellation to the centre
length_x = compact_tessellation{1,n}(8,1) - compact_tessellation{1,1}(9,1);
length_y = compact_tessellation{m,1}(14,2) - compact_tessellation{1,1}(9,2);
for i = 1:m
    for j = 1:n
        t = [n/2-length_x/2,m/2-length_y/2];
        compact_tessellation{i,j} = transform_square(compact_tessellation{i,j}(1:16,:),0,t,compact_tessellation{i,j}(9,:));
    end
end
end

