function edge_diff = edge_diff(tessellation_transformed,m,n)
% Calculate the edge length difference bewteen two coresponding unit.
tessellation = {};
for i = 1:m
    for j = 1:n
        tessellation{i,j} = [tessellation_transformed((i-1)*n*16+16*(j-1)+1:(i-1)*n*16+16*j,:)];
    end
end
index = [1 2;2 3;3 4;4 1;5 6;6 7;7 8;8 5;9 10;10 11;11 12;12 9;13 14;14 15;15 16;16 13];
length_total = 0;
for i = 1:m-1
    for j = 1:n-1
        for l = 1:16
            nodes1 = tessellation{i,j}([index(l,:)],:);
            nodes2 = tessellation{i+1,j}([index(l,:)],:);
            nodes3 = tessellation{i,j+1}([index(l,:)],:);
            length1 = length_calculate(nodes1);
            length2 = length_calculate(nodes2);
            length3 = length_calculate(nodes3);
            length_total = length_total + (length1 - length2)^2 + (length1 - length3)^2;
        end
    end
end
for l = 1:16
    nodes1 = tessellation{m,n}([index(l,:)],:);
    nodes2 = tessellation{m-1,n}([index(l,:)],:);
    nodes3 = tessellation{m,n-1}([index(l,:)],:);
    length1 = length_calculate(nodes1);
    length2 = length_calculate(nodes2);
    length3 = length_calculate(nodes3);
    length_total = length_total + (length1 - length2)^2 + (length1 - length3)^2;
end
edge_diff = length_total;
end