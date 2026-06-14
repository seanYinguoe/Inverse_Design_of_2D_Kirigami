function tessellation_optimized = tessellation_optimization(tessellation_transformed,s,r,p)
% using fmincon to optimizee the configuration
%% Define the objective function
l = 1;
[m , n] = size(tessellation_transformed);
M = (m-1)*(n-1)*2 + (m-1)+(n-1);
%fun = @(tessellation_optimized) 1/M * (l/pi * (angle_diff(tessellation_optimized,m,n)) + edge_diff(tessellation_optimized,m,n)); 
fun = @(tessellation_optimized) 1/M * (l/pi * (angle_diff(tessellation_optimized,m,n)) + edge_diff(tessellation_optimized,m,n)); 
%% Define the initial configuration
tessellation_initial = zeros(m*n*16,2);
for i = 1:m
    for j = 1:n
        tessellation_initial((i-1)*n*16+1+(j-1)*16:j*16+(i-1)*n*16,:) = tessellation_transformed{i,j};
    end
end
%% define the nonlinear condition
%nonlcon = @(tessellation_optimized) tessellation_con(tessellation_optimized,s,m,n,r);
if p == 1
    nonlcon = @(tessellation_optimized) rigid(tessellation_optimized,s,m,n,r);
elseif p == 2
        nonlcon = @(tessellation_optimized) nonrigid(tessellation_optimized,s,m,n,r);
end
A = []; % linear inequality constraints
b = []; % linear inequality constraintsyo
Aeq = zeros(2*m*n*16,2*m*n*16); % linear equality constraints
beq = zeros(2*m*n*16,1); % linear equality constraints
lb = []; % Lower bounds
ub = []; % Upper bounds
% connection of nodes
index1 = [4 10 6 16 8 3 15 2];
index2 = [7 13 11 1 9 14 12 5];
l = 1;
for i = 1:m
    for j = 1:n
        for k = 1:4  
            Aeq(l,index1(k)+(i-1)*n*16+16*(j-1)) = 1;  % x cooridinates
            Aeq(l,index2(k)+(i-1)*n*16+16*(j-1)) = -1;
            Aeq(l+1,index1(k)+(i-1)*n*16+16*(j-1)+m*n*16) = 1; % y coordinates
            Aeq(l+1,index2(k)+(i-1)*n*16+16*(j-1)+m*n*16) = -1;
            l = l+2;
        end
    end
end
for i = 1:m
    for j = 1:n-1
        for k = 5:6
            Aeq(l,index1(k)+(i-1)*n*16+16*(j-1)) = 1;
            Aeq(l,index2(k)+(i-1)*n*16+16*j) = -1;
            Aeq(l+1,index1(k)+(i-1)*n*16+16*(j-1)+m*n*16) = 1;
            Aeq(l+1,index2(k)+(i-1)*n*16+16*j+m*n*16) = -1;
            l = l+2;
        end
    end
end
for i = 1:m-1
    for j = 1:n
        for k = 7:8
            Aeq(l,index1(k)+(i-1)*n*16+16*(j-1)) = 1;
            Aeq(l,index2(k)+i*n*16+16*(j-1)) = -1;
            Aeq(l+1,index1(k)+(i-1)*n*16+16*(j-1)+m*n*16) = 1;
            Aeq(l+1,index2(k)+i*n*16+16*(j-1)+m*n*16) = -1;
            l = l+2;
        end
    end
end
% Define the options for the solver
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point','MaxFunEvals',8000);
%options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point',  'StepTolerance', 1e-6);
% Call the fmincon function to minimize the objective function
[tessellation_optimized1, fval] = fmincon(fun, tessellation_initial, A, b, Aeq, beq, lb, ub, nonlcon, options);
tessellation_optimized = cell(m,n);
for i = 1:m
    for j = 1:n
        tessellation_optimized{i,j} = [tessellation_optimized1((i-1)*n*16+16*(j-1)+1:(i-1)*n*16+16*j,:)];
    end
end
