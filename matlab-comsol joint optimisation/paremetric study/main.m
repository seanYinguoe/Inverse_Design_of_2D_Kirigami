%function bestchrom=GGA
%% reset programme

%clc
%clear
warning off      
format short g   
% study the parametric sensitivity of kirigami tessellation
% the width of cuts
% the gap between the cuts

%% build the geometry model
m = 2;
n = 2;          % m:the number of units in y direction, n:the number of unit in x direction
% create the intial tessellation
%temp = tessellation_intial(1:m/2,1:n/2); % create the initial tessellation
temp = tessellation_deployment(m,n,1,0);
%temp = temp(1:m/2,1:n/2);
for i = 1:m/2
    for j = 1:n/2
        temp{i,j} = temp{i,j}([1,2,3,4,5,8,9,13,14],:);
    end
end
nodes = derive_nodes(temp);
tessellation = zeros(m*n*9/4,2);
for i = 1:m/2
    for j = 1:n/2
        tessellation((i-1)*n/2*9+1+(j-1)*9:(i-1)*n/2*9+(j-1)*9+9,:) = nodes{i,j};
    end
end
temp = [tessellation(:,1);tessellation(:,2)];  % initial cut pattern
% derive the optimised variables
index = [];
% x,y coordinates
for i = 1:m/2
    for j = 1:n/2
        index = horzcat(index,[1+(i-1)*n/2*9+(j-1)*9]);
        index = horzcat(index,[1+(i-1)*n/2*9+(j-1)*9+m*n/4*9]);
    end
end
if n > 2
    for i = 1:m/2
        for j = 1:(n/2-1)
            index = horzcat(index,[4+(i-1)*n/2*9+(j-1)*9]);
            index = horzcat(index,[4+(i-1)*n/2*9+(j-1)*9+m*n/4*9]);
        end
    end
end
if m > 2
    for i = 1:(m/2-1)
        for j = 1:n/2
            index = horzcat(index,[2+(i-1)*n/2*9+(j-1)*9]);
            index = horzcat(index,[2+(i-1)*n/2*9+(j-1)*9+m*n/4*9]);
        end
    end
end
if m > 2 && n > 2
    for i = 1:m/2-1
        for j = 1:n/2-1
            index = horzcat(index,[3+(i-1)*n/2*9+(j-1)*m/2*9]);
            index = horzcat(index,[3+(i-1)*n/2*9+(j-1)*m/2*9+m*n/4*9]);
        end
    end
end
% x coordinate
for j = 1:(n/2)
    index = horzcat(index,[5+(j-1)*9]);
    index = horzcat(index,[2+(j-1)*9+(n/2)*9*(m/2-1)]);
end
if n > 2
    for j = 1:(n/2-1)
        index = horzcat(index,[6+(j-1)*9]);
        index = horzcat(index,[3+(j-1)*9+(n/2)*9*(m/2-1)]);
    end
end
% y coordinate
for i = 1:(m/2)
    index = horzcat(index,[8+(i-1)*9*n/2+m*n/4*9]);
    index = horzcat(index,[4+(i-1)*9*n/2+9*(n/2-1)+m*n/4*9]);
end
if m > 2
    for i = 1:(m/2-1)
        index = horzcat(index,[9+(i-1)*9+m*n/4*9]);
        index = horzcat(index,[3+(i-1)*9+9*(n/2-1)+m*n/4*9]);
    end
end
x0 = temp(index');
x_initial = x0;
 
%% set the material properties(buckling-induced Kirigami)

density = 2700;                 % density of material
E = 4.33e9;                       % the elastic module of material
nu = 0.4;                       % possion ratio
properties = [density,E,nu];    % set the property of material

%% set the parameters

displacement = 0.3;                 % set the displacement
g = 0.03;                             % set the gap between cuts
w = 0.03;                             % set the width of cuts

%% get the results
res = parametric(m,n,x0,displacement,properties,g,w); 








