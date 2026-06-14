function tessellation_transformed = tessellation_deployment(m,n,length,Rotate_angels)
% Deployment of squares tessellation
% intialize the tessellation
% num = 4;  % the number of tessellation units in each axis 
% length = 1;  % the length of each units
% Rotate_angels = pi/4; % the rotation angles of squares tessellation
tessellation = create_tessellation(length,m,n); % create the squares tessellation
% Rotate each square in a single unit
[m,n] = size(tessellation);
tessellation_transformed = cell(m,n);
for i = 1:m
    for j = 1:n
        t1 = [length/2*(sin(Rotate_angels)+cos(Rotate_angels)-1)*(j-1)*2,length/2*(sin(Rotate_angels)+cos(Rotate_angels)-1)*(i-1)*2];
        t = [length/2*(1-cos(Rotate_angels)),0]; % the traslation vector
        square1 = tessellation{i,j}([1,2,3,4],:);
        square2 = tessellation{i,j}([5,6,7,8],:);
        square3 = tessellation{i,j}([9,10,11,12],:);
        square4 = tessellation{i,j}([13,14,15,16],:);

        square_transformed1 = transform_square(square1, -Rotate_angels, -t+t1,tessellation{i,j}(4,:));
        square_transformed2 = transform_square(square2, Rotate_angels, -t+t1,tessellation{i,j}(7,:));
        square_transformed3 = transform_square(square3, -Rotate_angels, t+t1,tessellation{i,j}(10,:));
        square_transformed4 = transform_square(square4, Rotate_angels, t+t1,tessellation{i,j}(13,:));

        tessellation_transformed{i,j} = [square_transformed1([1:4],:);
            square_transformed2([1:4],:);
            square_transformed3([1:4],:);
            square_transformed4([1:4],:)];
    end
end
% move tessellation to centre (0,0)
xc = (tessellation_transformed{1,1}(9,1) + tessellation_transformed{1,n}(8,1))/2;
yc = (tessellation_transformed{1,1}(12,2) + tessellation_transformed{m,1}(15,2))/2;
for i = 1:m
    for j = 1:n
        %tessellation_transformed{i,j} = [tessellation_transformed{i,j}(:,1) - xc,tessellation_transformed{i,j}(:,2) - yc];
        tessellation_transformed{i,j} = [tessellation_transformed{i,j}(:,1) - xc,tessellation_transformed{i,j}(:,2) - yc];
    end
end
end
% Plot deployed tessellation
% figure(1);
% plot_tessellation(tessellation_transformed);