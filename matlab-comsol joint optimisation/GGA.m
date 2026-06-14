%function bestchrom=GGA
%% reset programme
%clc
%clear
warning off      
format short g   

%% parameters for GA
m = 4;
n = 4;          % m:the number of units in y direction, n:the number of unit in x direction
maxgen = 6;                                      % maximum evolutionary 
r = 3;                                         % a of ellips
sizepop=15;                                    % size of population
pcross = [0.6];                                % crossover rate
pmutation = [0.01];                            % mutation rate
num_v = 2*(n-1)*(m-1) + (n-1)*2 + (m-1)*2;        % the number of variables
lenchrom=ones(num_v,1);                           % length of variables
% create the intial tessellation
temp = tessellation_compacted(1:m/2,1:n/2); % create the initial tessellation
%temp = tessellation_deployment(m,n,1,0);
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
% leftover coordinate
% index = horzcat(index,[1+(n/2-1)*9]);
% index = horzcat(index,[1+(n/2-1)*9+m*n/4*9]);
% index = horzcat(index,[2+(n/2-1)*9]);
% index = horzcat(index,[5+(n/2-1)*9]);
%x0 = temp([1:6,10:11,14,19,22,26,28,31]);
x0 = temp(index');
x_initial = x0;
bound = [x0-0.1 x0+0.1];                 % bound for varaibles
% bound([7,8,9],:) = -2;
% bound([12,13,15],:) = 0;
% bound([23,24,25,32,33,34],:) = -1;
% bound([20,21,27,29,30,36],:) = 0;       
%% set the material properties
density = 2700;                 % density of material
E = 70e9;                       % the elastic module of material
nu = 0.2;                       % possion ratio
properties = [density,E,nu];    % set the property of material

%% set the load condition
d = 0.1;                 % set the displacement

%% initial population
individuals=struct('shape',[],'fitness',zeros(1,sizepop),'chrom',[]);   % structure=fitness+chrom
avgfitness=[];                                               % average fitness in population
bestfitness=[];                                              % bestfitness in population
bestchrom=[];                                                % bestchrome in population
% initial population 
for i=1:sizepop
    if i == 1
        individuals.chrom(i,:) = x0;                   % put initial cut pattern to intial population
    else
        individuals.chrom(i,:) = Code(lenchrom,bound);   % produce random population
    end                                                
    x = individuals.chrom(i,:);
    [individuals.fitness(i) individuals.shape(i,:)] = livelink(m,n,x,d,properties);    % calculate the individual fitness
end

    
% find the best chrome
[bestfitness bestindex] = min(individuals.fitness);      % find the best fitness and corresponding chrome
bestchrom = individuals.chrom(bestindex,:);              % best chrome
bestshape = individuals.shape(bestindex,:);              % best shape
avgfitness = sum(individuals.fitness)/sizepop;           % average fitness
% record best fitness and average fitness in every iterations
trace_fitness = [];
% record best boundary shape in every iterations
trace_shape = [];
% record best chrom in every iterations
trace_chrom = [];

%% start
for i=1:maxgen
    disp('start');
    % select
    individuals=Select(individuals,sizepop);
    avgfitness=sum(individuals.fitness)/sizepop;
    % crossover
    individuals.chrom=Cross(pcross,lenchrom,individuals.chrom,sizepop,bound);
    % mutation
    individuals.chrom=Mutation(pmutation,lenchrom,individuals.chrom,sizepop,[i maxgen],bound);
    % calculate fitness
    for j=1:sizepop
        disp('calculate times');
        disp([i j]);
        x=individuals.chrom(j,:);
        [individuals.fitness(j) individuals.shape(j,:)]=livelink(m,n,x,d,properties); 
    end
    
    % find the best and worst fitness and corrsponding chrome
    [newbestfitness,newbestindex]=min(individuals.fitness);
    [worestfitness,worestindex]=max(individuals.fitness);
    % use new bestfitness replace the last iterations
    if bestfitness>newbestfitness
        bestfitness=newbestfitness;
        bestchrom=individuals.chrom(newbestindex,:);
        bestshape=individuals.shape(newbestindex,:);
    end
    individuals.chrom(worestindex,:)=bestchrom;
    individuals.fitness(worestindex)=bestfitness;
    individuals.shape(worestindex,:)=bestshape;

    avgfitness=sum(individuals.fitness)/sizepop;
    trace_fitness=[trace_fitness;avgfitness bestfitness];
    trace_chrom=[trace_chrom;bestchrom];
    trace_shape=[trace_shape;bestshape];
    disp('end');
    disp(' ');
    disp(' ');
end
save trace_fitness trace_fitness
save trace_shape trace_shape

%% result
figure(1)
[r c]=size(trace_fitness);
plot([1:r]',trace_fitness(:,1),'r-',[1:r]',trace_fitness(:,2),'b--');  
title(['cost function  ' 'iterations' num2str(maxgen)],'fontsize',12);
xlabel('iteration number ','fontsize',12);ylabel('cost function','fontsize',512);
legend('average value','best value','fontsize',12);
ylim()   
grid on
 
x=bestchrom;
disp('best value');
figure(2)
livelink(m,n,x,d,properties);

%% plot figure
% plot the boundary shape in iterations
[fitness_initial shape_initial] = livelink(m,n,x_initial,d,properties);
figure(3)
% initial boundary shape
[~,num] = size(shape_initial);
x0 = shape_initial(1:num/2);
y0 = shape_initial(num/2+1:num);
p = polyfit(x0,y0,9); % fit the a ploynomial function
x = linspace(-n/2-d,n/2+d);
y = polyval(p,x);
initial_shape = plot(x,y,'DisplayName',['initial shape'],'Color',[0 0.4470 0.7410],'LineWidth', 3,'LineStyle','-.');
initial_shape.Color(4) = 0.3;  % change the transparency of the curve
% iterative boundary shape
for i = 1:maxgen
    hold on
    x0 = trace_shape(i,1:num/2);
    y0 = trace_shape(i,num/2+1:num);
    p = polyfit(x0,y0,3); % fit the a ploynomial function
    x = linspace(-n/2-d,n/2+d);
    y = polyval(p,x);
    expr = ['shape_' num2str(i)];
    expr = plot(x,y,'DisplayName',['iteration ' num2str(i)],'Color', [0 0.4470 0.7410],'LineWidth', 3);
    expr.Color(4) = 0.1*i;  % change the transparency of the curve
    xlim([-n/2-0.5,n/2+0.5]);
    ylim([-m/2-0.5,-m/2+0.5]);
end
% target boundary shape
gap = shape_initial(1,num/2+1)+1/2*sqrt(r.^2-(-n/2-d).^2);
x_target = linspace(-n/2-d,n/2+d);
y_target = -1/sqrt(6)*sqrt(r.^2-x_target.^2)+gap;
target_shape = plot(x_target,y_target,'DisplayName',['target shape'],'Color',[0.6350 0.0780 0.1840],'LineWidth', 3,'LineStyle','-.');
legend();
hold off
figure(4)

