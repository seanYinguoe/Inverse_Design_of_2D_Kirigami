function [C1,C2] = cuts_rigid(nodes,d)
%CUTS_RIGID Original cut paths, shortened by d at the ligament ends.
[m,n]=size(nodes);

k1 = 1;
k2 = 1;
C1 = [];
C2 = [];
% build the whole cuts
for i = 1:m
    for j = 1:n
        temp = nodes{i,j}(9,:) - nodes{i,j}(4,:);
        alpha = d/norm(temp);
        C1(k1,:) = [nodes{i,j}(9,:)-alpha*temp,nodes{i,j}(4,:)+alpha*temp];
        k1 = k1+1;
    end
end
for i = 1:m-1
    for j = 1:n-1
        temp = nodes{i,j}(2,:) - nodes{i,j+1}(2,:);
        alpha = d/norm(temp);
        C1(k1,:) = [nodes{i,j}(2,:)-alpha*temp,nodes{i,j+1}(2,:)+alpha*temp];
        k1 = k1+1;
    end
end
for i = 1:m-1
    for j = 1:n
        temp = nodes{i,j}(1,:) - nodes{i+1,j}(6,:);
        alpha = d/norm(temp);
        C1(k1,:) = [nodes{i,j}(1,:)-alpha*temp,nodes{i+1,j}(6,:)+alpha*temp];
        k1 = k1+1;
    end
end
for i = 1:m
    for j = 1:n-1
        temp = nodes{i,j}(3,:) - nodes{i,j}(7,:);
        alpha = d/norm(temp);
        C1(k1,:) = [nodes{i,j}(3,:)-alpha*temp,nodes{i,j}(7,:)+alpha*temp];
        k1 = k1+1;
    end
end
% build the half cut
for j = 1:n
    temp = nodes{1,j}(6,:) - nodes{1,j}(5,:);
    alpha = d/norm(temp);
    C2(k2,:) = [nodes{1,j}(6,:) - alpha*temp,nodes{1,j}(5,:) - alpha*temp];
    k2 = k2+1;
end
for j = 1:n
    temp = nodes{m,j}(2,:) - nodes{m,j}(1,:);
    alpha = d/norm(temp);
    C2(k2,:) = [nodes{m,j}(2,:) + alpha*temp,nodes{m,j}(1,:) + alpha*temp];
    k2 = k2+1;
end
for i = 1:m-1
    temp = nodes{i,1}(10,:) - nodes{i,1}(2,:);
    alpha = d/norm(temp);
    C2(k2,:) = [nodes{i,1}(10,:) + alpha*temp,nodes{i,1}(2,:) + alpha*temp];
    k2 = k2+1;
end
for i = 1:m-1
    temp = nodes{i,n}(3,:) - nodes{i,n}(2,:);
    alpha = d/norm(temp);
    C2(k2,:) = [nodes{i,n}(3,:) + alpha*temp,nodes{i,n}(2,:) + alpha*temp];
    k2 = k2+1;
end

end
