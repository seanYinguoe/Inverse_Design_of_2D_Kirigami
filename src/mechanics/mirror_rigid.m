function nodes_s = mirror_rigid(tessellation)
%% duplicate the effective nodes for build the kirigami pattern
% get the full kirigami tesselltion from 1/4 nodes
[m n] = size(tessellation);
nodes_s = tessellation;
for i = 1:m
    for j = 1:n
        if i == m
            nodes_s{m,j}([2,3,10],2) = 0;
        end
        if j == n
            nodes_s{i,n}([3,4,7],1) = 0;
        end
        nodes_s{2*m+1-i,j} = [nodes_s{i,j}([6,5,7,4,2,1,3,10,9,8],1) , -nodes_s{i,j}([6,5,7,4,2,1,3,10,9,8],2)]; % duplicate symmetric nodes about x axis
        nodes_s{i,2*n+1-j} = [-nodes_s{i,j}([1,2,10,9,5,6,8,7,4,3],1),nodes_s{i,j}([1,2,10,9,5,6,8,7,4,3],2)]; % duplicate symmetric nodes about y axis
        nodes_s{2*m+1-i,2*n+1-j} = [-nodes_s{i,j}([6,5,8,9,2,1,10,3,4,7],1),-nodes_s{i,j}([6,5,8,9,2,1,10,3,4,7],2)]; % duplicate symmetric nodes about zero point
    end
end

end
