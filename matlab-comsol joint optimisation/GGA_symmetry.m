%function bestchrom=GGA
%% reset programme
warning off      
format short g  

%% parameters for GA
m = 2;
n = 4;          % m:the number of units in y direction, n:the number of unit in x direction
maxgen = 1;     % maximum evolutionary 
gen = 0;        % evolutionary
s = 2;          % s = 1,ellispe;s = 2,vase;s = 3,wavy
%r = 3.2*sqrt(1);                                         % a of ellips
r = 3.2*sqrt(1);
sizepop = 1;                                    % size of population
pcross = [0.6];                                % crossover rate
pmutation = [0.01];                            % mutation rate
num_v = 2*(n-1)*(m-1) + (n-1)*2 + (m-1)*2;        % the number of variables
lenchrom=ones(num_v,1);                           % length of variables

%% create the intial tessellation
temp = tessellation_compacted(1:m/2,1:n/2); % create the initial tessellation based on the kinematic optimised pattern
%temp = tessellation_intial(1:m/2,1:n/2);
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
temp = real(temp);

%% derive the optimised variables
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
            %index = horzcat(index,[9+(i-1)*n/2*9+(j-1)*9]);
            %index = horzcat(index,[9+(i-1)*n/2*9+(j-1)*9+m*n/4*9]);
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

% other nodes
x0 = temp(index');
x_initial = x0;

%% define the boundary condition 
bound = [x0-0.08 x0+0.08];                 % bound for varaibles  

%% define the material properties
density = 2700;                 % density of material
E = 70e9;                       % the elastic module of material
nu = 0.2;                       % possion ratio
properties = [density,E,nu];    % set the property of material

%% set the load/displacement condition
d = 0.5;                 % set the displacement

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
    try
        [individuals.fitness(i) individuals.shape(i,:)] = livelink(m,n,x,d,properties,s);    % calculate the individual fitness
    catch
        x = x0;
        [individuals.fitness(i) individuals.shape(i,:)] = livelink(m,n,x,d,properties,s);
    end
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
% set the error thresholds
max_count = 4;
stable_count = 0;
error_threshold = 0.005;
flag = 1;
i=1;
%% start GA
while flag
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
        try
            [individuals.fitness(j) individuals.shape(j,:)]=livelink(m,n,x,d,properties,s); 
        catch
            individuals.chrom(j,:)=bestchrom;
            x=individuals.chrom(j,:);
            [individuals.fitness(j) individuals.shape(j,:)]=livelink(m,n,x,d,properties,s); 
        end
    end
    
    % find the best and worst fitness and corrsponding chrome
    [newbestfitness,newbestindex]=min(individuals.fitness);
    [worestfitness,worestindex]=max(individuals.fitness);
    % use new bestfitness replace the last iterations
    if bestfitness>newbestfitness
        bestfitness=newbestfitness;
        bestchrom=individuals.chrom(newbestindex,:);
        bestshape=individuals.shape(newbestindex,:);
    else
        stable_count = stable_count + 1;
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
    %| stable_count >= max_count
    if (bestfitness < error_threshold) | (i >=maxgen)
        flag = 0;
    end
    i = i+1;
    gen = gen+1;
end
save trace_fitness trace_fitness
save trace_shape trace_shape

%% result
figure(1)
[row col]=size(trace_fitness);
plot([1:row]',trace_fitness(:,1),'r-',[1:row]',trace_fitness(:,2),'b--');  
title(['cost function  '],'fontsize',12);
xlabel('iteration number ','fontsize',12);ylabel('cost function','fontsize',12);
legend('average value','best value','fontsize',12);
ylim()   
grid off
 
x=bestchrom;
disp('best value');
figure(2)
livelink(m,n,x,d,properties,s);

%% plot figure
% plot the boundary shape in iterations
[fitness_initial shape_initial] = livelink(m,n,x_initial,d,properties,s);
figure(3)
xlim([-3 3]);
ylim([-4 -1.5]);

%% ellipse initial boundary shape
figure(4)
xlabel('X(cm)', 'FontSize', 16) 
ylabel('Y(cm)', 'FontSize', 16)
title('Boundary shape', 'FontSize', 16) 
set(gca, 'FontSize', 14)


if s == 1
    [~,num] = size(shape_initial);
    x0 = shape_initial(1:num/2);
    y0 = shape_initial(num/2+1:num);
    p = polyfit(x0,y0,7); % fit the a ploynomial function
    syms x
    p1 = poly2sym(p,x);
    gap = subs(p1,x,-2.5)+1/sqrt(4)*sqrt(1*r.^2-(-n/2-d).^2);
    x = linspace(-n/2-d,n/2+d);
    y = polyval(p,x) - gap;
    % plot initial boundary shape
    % initial_points = scatter(x0,y0-gap,'filled','Color',[0 0.4470 0.7410],'MarkerFaceAlpha', 0.3);
    hold on
    initial_shape = plot(x,y,'DisplayName',['initial shape'],'Color',[0 0.4470 0.7410],'LineWidth', 3,'LineStyle','-.');
    initial_shape.Color(4) = 0.3;  % change the transparency of the curve
    % iterative boundary shape
    for i = 1:5:maxgen
        hold on
        x0 = trace_shape(i,1:num/2);
        y0 = trace_shape(i,num/2+1:num);
        p = polyfit(x0,y0,7); % fit the a ploynmial function
        x = linspace(-n/2-d,n/2+d);
        syms x1
        p1 = poly2sym(p,x1);
        gap = subs(p1,x1,-2.5)+1/sqrt(4)*sqrt(1*r.^2-(-n/2-d).^2);
        y = polyval(p,x) - gap;
        expr = ['shape_' num2str(i)];
        expr = plot(x,y,'DisplayName',['iter. ' num2str(i)],'Color', [0 0.4470 0.7410],'LineWidth', 3);
        expr.Color(4) = 1/maxgen*i;  % change the transparency of the curve
    end
    % target boundary shape
    x_target = linspace(-n/2-d,n/2+d);
    y_target = -1/sqrt(4)*sqrt(1*r.^2-x_target.^2);
    target_shape = plot(x_target,y_target,'DisplayName',['target shape'],'Color',[0.6350 0.0780 0.1840],'LineWidth', 3,'LineStyle','-.');
    legend();
    hold off

    %% vase initial boundary shape
elseif s == 2
    [~,num] = size(shape_initial);
    xi = shape_initial(1:num/2);
    yi = shape_initial(num/2+1:num);
    p = polyfit(xi,yi,7); % fit the a ploynomial function
    gap = polyval(p,-n/2-0.5) - (sqrt(1/3*r^2-1/3*(-n/2-0.5)^2)-r);
    x = linspace(-n/2-d,n/2+d);
    y = polyval(p,x) - gap;
    initial_shape = plot(x,y,'DisplayName',['initial shape'],'Color',[0 0.4470 0.7410],'LineWidth', 3,'LineStyle','-.');
    initial_shape.Color(4) = 0.3;  % change the transparency of the curve
   % iterative boundary shape
    for i = 1:5:maxgen
        hold on
        x0 = trace_shape(i,1:num/2);
        y0 = trace_shape(i,num/2+1:num);
        p = polyfit(x0,y0,7); % fit the a ploynmial function
        x = linspace(-n/2-d,n/2+d);
        y = polyval(p,x);
        gap = polyval(p,-n/2-d) - (sqrt(1/3*r^2-1/3*(-n/2-d)^2)-r);
        y = polyval(p,x) - gap;
        expr = ['shape_' num2str(i)];
        expr = plot(x,y,'DisplayName',['iteration ' num2str(i)],'Color', [0 0.4470 0.7410],'LineWidth', 3);
        expr.Color(4) = 1/maxgen*i;  % change the transparency of the curve
    end
    % target boundary shape
    x_target = linspace(-n/2-d,n/2+d);
    y_target = sqrt(1/3*r.^2-1/3*(x_target).^2)-r;
    target_shape = plot(x_target,y_target,'DisplayName',['target shape'],'Color',[0.6350 0.0780 0.1840],'LineWidth', 3,'LineStyle','-.');
    legend();
    hold off

elseif s == 3
    % sinwave initial boundary shape
    [~,num] = size(shape_initial);
    x0 = shape_initial(1:num/2);
    y0 = shape_initial(num/2+1:num);
    p = polyfit(x0,y0,7); % fit the a ploynomial function
    gap = shape_initial(1,num/2+1)+(0.2*cos(0.8*pi*(-n/2-0.5-1/0.8)) + 1.5);
    x = linspace(-n/2-d,n/2+d);
    y = polyval(p,x) - gap;
    initial_shape = plot(x,y,'DisplayName',['initial shape'],'Color',[0 0.4470 0.7410],'LineWidth', 3,'LineStyle','-.');
    initial_shape.Color(4) = 0.3;  % change the transparency of the curve
    % iterative boundary shape
    for i = 1:5:gen
        hold on
        x0 = trace_shape(i,1:num/2);
        y0 = trace_shape(i,num/2+1:num);
        p = polyfit(x0,y0,7); % fit the a ploynmial function
        x = linspace(-n/2-d,n/2+d);
        y = polyval(p,x);
        gap = polyval(p,-n/2-d)+(0.2*cos(0.8*pi*(-n/2-0.5-1/0.8)) + 1.5);
        y = polyval(p,x) - gap;
        expr = ['shape_' num2str(i)];
        expr = plot(x,y,'DisplayName',['iteration ' num2str(i)],'Color', [0 0.4470 0.7410],'LineWidth', 3);
        expr.Color(4) = 1/gen*i;  % change the transparency of the curve
    end
    % target boundary shape
    x_target = linspace(-n/2-d,n/2+d);
    y_target = -(0.2*cos(0.8*pi*(x_target-1/0.8)) + 1.5);
    target_shape = plot(x_target,y_target,'DisplayName',['target shape'],'Color',[0.6350 0.0780 0.1840],'LineWidth', 3,'LineStyle','-.');
    legend();
    hold off
end





