function nodes_s = derive_nodes(tessellation)
%% duplicate the effective nodes for build the kirigami pattern
% get the 1/4 nodes symmetry
[m n] = size(tessellation);
nodes_s = tessellation;
for i = 1:m
    for j = 1:n
        if i == m
            nodes_s{m,j}([2,3,9],2) = 0;
        end
        if j == n
            nodes_s{i,n}([3,4,6],1) = 0;
        end
        nodes_s{2*m+1-i,j} = [nodes_s{i,j}([1,5,6,4,2,3,9,8,7],1) , -nodes_s{i,j}([1,5,6,4,2,3,9,8,7],2)]; % duplicate symmetric nodes about x axis
        nodes_s{i,2*n+1-j} = [-nodes_s{i,j}([1,2,9,8,5,7,6,4,3],1),nodes_s{i,j}([1,2,9,8,5,7,6,4,3],2)]; % duplicate symmetric nodes about y axis
        nodes_s{2*m+1-i,2*n+1-j} = [-nodes_s{i,j}([1,5,7,8,2,9,3,4,6],1),-nodes_s{i,j}([1,5,7,8,2,9,3,4,6],2)]; % duplicate symmetric nodes about zero point
    end
end

