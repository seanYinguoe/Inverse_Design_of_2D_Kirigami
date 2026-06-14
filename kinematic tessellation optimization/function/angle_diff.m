function angle_diff = angle_diff(tessellation_transformed,m,n)
% Calculate the angle difference bewteen two coresponding unit.
tessellation = {};
for i = 1:m
    for j = 1:n
        tessellation{i,j} = [tessellation_transformed((i-1)*n*16+16*(j-1)+1:(i-1)*n*16+16*j,:)];
    end
end
angle_total = 0;
index = [1 2 3;2 3 4;3 4 1;4 1 2;5 6 7;6 7 8;7 8 5;8 5 6;9 10 11;10 11 12;
    11 12 9;12 9 10;13 14 15;14 15 16;15 16 13;16 13 14];
for i = 1:m-1
    for j = 1:n-1
        for l = 1:16
            nodes1 = tessellation{i,j}([index(l,:)],:);
            nodes2 = tessellation{i+1,j}([index(l,:)],:);
            nodes3 = tessellation{i,j+1}([index(l,:)],:);
            angle_calculate1 = (angle_calculate(nodes1) - angle_calculate(nodes2))^2;
            angle_calculate2 = (angle_calculate(nodes1) - angle_calculate(nodes3))^2;
            angle_total = angle_total + angle_calculate1 + angle_calculate2;
        end
    end
end
for l = 1:16
    nodes1 = tessellation{m,n}([index(l,:)],:);
    nodes2 = tessellation{m-1,n}([index(l,:)],:);
    nodes3 = tessellation{m,n-1}([index(l,:)],:);
    angle_calculate1 = (angle_calculate(nodes1) - angle_calculate(nodes2))^2;
    angle_calculate2 = (angle_calculate(nodes1) - angle_calculate(nodes3))^2;
    angle_total = angle_total + angle_calculate1 + angle_calculate2;
end
angle_diff = angle_total;
end