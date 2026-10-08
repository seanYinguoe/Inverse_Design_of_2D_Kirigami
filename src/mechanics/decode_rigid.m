function nodes = decode_rigid(tessellation,m,n)
%DECODE_RIGID Original straight-cut parameterisation for a symmetric 2-by-4 grid.
assert(m==2&&n==4&&numel(tessellation)==12, ...
    'The rigid mechanical parameterisation supports a symmetric 2-by-4 grid only.');
temp = cell(m/2,n/2);

% derive the optimised variables
index = [];

% x coordinate
for j = 1:(n/2)
    index = horzcat(index,[5+(j-1)*10]);
    index = horzcat(index,[2+(j-1)*10+(n/2)*10*(m/2-1)]);
end
if n > 2
    for j = 1:(n/2-1)
        index = horzcat(index,[7+(j-1)*10]);
        index = horzcat(index,[3+(j-1)*10+(n/2)*10*(m/2-1)]);
    end
end

% y coordinate
for i = 1:(m/2)
    index = horzcat(index,[9+(i-1)*10*n/2+m*n/4*10]);
    index = horzcat(index,[4+(i-1)*10*n/2+10*(n/2-1)+m*n/4*10]);
end
if m > 2
    for i = 1:(m/2-1)
        index = horzcat(index,[10+(i-1)*10+m*n/4*10]);
        index = horzcat(index,[3+(i-1)*10+10*(n/2-1)+m*n/4*10]);
    end
end
% nodes extension from num_v variables to m*n/4*10*2 variables
x = zeros(m*n*20/4,1);
x(index) = tessellation(1:8);
for i = 1:(m/2)
    x([8,9,10]+(n/2)*(i-1)*10) = -n/2;
end
for i = 1:m/2
    for j = 1:(n/2-1)
        x([18 20]+10*(j-1)+(i-1)*(n/2)*10) = x([7 3]+10*(j-1)+(i-1)*(n/2)*10);
        x([18 20]+10*(j-1)+10*m*n/4+(i-1)*(n/2)*10) = x([7 3]+10*(j-1)+10*m*n/4+(i-1)*(n/2)*10);
    end
end
for j = 1:n/2
    x([5 7 8]+10*(j-1)+10*m*n/4) = -m/2;
end
% add nodes 4, 9, 1, 6
node_4 = tessellation(9) * ([tessellation(6),0]-[tessellation(5),-m/2]) + [tessellation(5),-m/2];
node_9 = tessellation(10) * ([tessellation(6),0]-[tessellation(5),-m/2]) + [tessellation(5),-m/2];
node_1 = [tessellation(2) (tessellation(2)+n/2)/(node_4(1)+n/2)*(node_4(2)-tessellation(7))+tessellation(7)];
node_6 = tessellation(11) * (node_4 - [-n/2 tessellation(7)]) + [-n/2 tessellation(7)];
node_21 = [tessellation(4) (tessellation(4)-node_9(1))/(0-node_9(1))*(tessellation(8)-node_9(2))+node_9(2)];
node_26 = tessellation(12) * ([0 tessellation(8)] - node_9) + node_9;
x([1 21 4 24 19 39 6 26 16 36 11 31]) = [node_1 node_4 node_9 node_6 node_26 node_21];

% derive the nodes (m/2)*(n/2)
for i = 1:m/2
    for j = 1:n/2
        temp{i,j}(:,1) = [x((i-1)*n/2*10+10*(j-1)+1:(i-1)*n/2*10+10*j)];
        temp{i,j}(:,2) = [x(10*m*n/4+(i-1)*n/2*10+10*(j-1)+1:10*m*n/4+(i-1)*n/2*10+10*j)];
    end
end
nodes = mirror_rigid(temp); % dupicate nodes (m/2)*(n/2) to m*n

end
