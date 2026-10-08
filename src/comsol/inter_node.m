function nodes = inter_node(x1,x2,n,p)
% produce nodes bewteen two nodes with different intervals
% input: nodes cordinate x1 x2
%        number of nodes produced in between
%        clustering degree
nodes = [x1];
d = x2-x1;
x0 = x1;
for i = 1:n
    for j = 1:(1/p-1)
        nodes1 = x0 + j*d*(p)^i;
        nodes = [nodes nodes1];
    end
    x0 = x2 - d*p^i;
end
nodes = [nodes x2];
end
