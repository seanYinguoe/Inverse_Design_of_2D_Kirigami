function [res boundary_shape] = livelink_nonrigid(m,n,tessellation,displacement,properties)
% livelink matlab commands for comsol
%
% This model demonstrates how to use matlab optimize the kirigami pattern
% using the Livelink for matlab.
% units: m,
%% build geometry model in tessellation_optimization.m
temp = cell(m,n);
% derive the optimised variables
index = [];
% x,y coordinates
for i = 1:m
    for j = 1:n
        index = horzcat(index,[1+(i-1)*n*9+(j-1)*9]);
        index = horzcat(index,[1+(i-1)*n*9+(j-1)*9+m*n*9]);
    end
end
if n >= 2
    for i = 1:m
        for j = 1:n-1
            index = horzcat(index,[4+(i-1)*n*9+(j-1)*9]);
            index = horzcat(index,[4+(i-1)*n*9+(j-1)*9+m*n*9]);
        end
    end
end
if m >= 2
    for i = 1:m-1
        for j = 1:n
            index = horzcat(index,[2+(i-1)*n*9+(j-1)*9]);
            index = horzcat(index,[2+(i-1)*n*9+(j-1)*9+m*n*9]);
        end
    end
end
if m >= 2 && n >= 2
    for i = 1:m-1
        for j = 1:n-1
            index = horzcat(index,[3+(i-1)*n*9+(j-1)*9]);
            index = horzcat(index,[3+(i-1)*n*9+(j-1)*9+m*n*9]);
        end
    end
end
% x coordinate
for j = 1:n
    index = horzcat(index,[5+(j-1)*9]);
    index = horzcat(index,[2+(j-1)*9+n*9*(m-1)]);
end
if n >= 2
    for j = 1:n-1
        index = horzcat(index,[6+(j-1)*9]);
        index = horzcat(index,[3+(j-1)*9+n*9*(m-1)]);
    end
end
% y coordinate
for i = 1:m
    index = horzcat(index,[8+(i-1)*9*n+m*n*9]);
    index = horzcat(index,[4+(i-1)*9*n+9*(n-1)+m*n*9]);
end
if m >= 2
    for i = 1:m-1
        index = horzcat(index,[9+(i-1)*9*n+m*n*9]);
        index = horzcat(index,[3+(i-1)*9*n+9*(n-1)+m*n*9]);
    end
end
% nodes extension variables
x = zeros(m*n*18,1);
x(index) = tessellation;
for i = 1:m
    x([7,8,9]+n*(i-1)*9) = -n/2;
end
for i = 1:m
    x([3,4,6]+n*(i-1)*9+(n-1)*9) = n/2;
end
% the connected nodes
for i = 1:m
    for j = 1:n-1
        x([18 17 16]+9*(j-1)+(i-1)*(n)*9) = x([3 4 6]+9*(j-1)+(i-1)*(n)*9);
        x([18 17 16]+9*(j-1)+9*m*n+(i-1)*(n)*9) = x([3 4 6]+9*(j-1)+9*m*n+(i-1)*(n)*9);
    end
end
for i = 1:m-1
    for j = 1:n
        x([23 24 25]+9*(j-1)+(i-1)*(n)*9) = x([2 3 9]+9*(j-1)+(i-1)*(n)*9);
        x([23 24 25]+9*(j-1)+(i-1)*(n)*9+9*m*n) = x([2 3 9]+9*(j-1)+(i-1)*(n)*9+9*m*n);
    end
end
for j = 1:n
    x([5 6 7]+9*(j-1)+9*m*n) = -m/2;
end
for j = 1:n
    x([2 3 9]+9*(j-1)+9*m*n+9*(m-1)*n) = m/2;
end
for i = 1:m
    for j = 1:n
        temp{i,j}(:,1) = [x((i-1)*n*9+9*(j-1)+1:(i-1)*n*9+9*j)];
        temp{i,j}(:,2) = [x(9*m*n+(i-1)*n*9+9*(j-1)+1:9*m*n+(i-1)*n*9+9*j)];
    end
end
length_x = temp{1,n}(6,1) - temp{1,1}(7,1);
length_y = temp{m,1}(9,2) - temp{1,1}(7,2);
nodes = temp;
d = 0.03; %the gap between cuts
w = 0.02; %the width of cuts
f_nodes = [1 2 3 4];                % set the fixed nodes
%d_nodes = [177 178 179 180];        % set the displacement nodes 2*4
%d_nodes = [269 270 271 272];        %set the displacement nodes 2*6
% set the material properties
density = 2700;                 % density of material
E = 70e9;                       % the elastic module of material
nu = 0.2;                      % possion ratio
properties = [density,E,nu];    % set the property of material

% Define the geometry of cuts.
k1 = 1;
k2 = 1;
C1 = [];
C2 = [];
% build the whole cuts
for i = 1:m
    for j = 1:n
        temp1 = nodes{i,j}(8,:) - nodes{i,j}(1,:);
        temp2 = nodes{i,j}(4,:) - nodes{i,j}(1,:);
        alpha1 = d/norm(temp1);
        alpha2 = d/norm(temp2);
        C1(k1,:) = [nodes{i,j}(8,:)-2*alpha1*temp1,nodes{i,j}(1,:),nodes{i,j}(4,:)-2*alpha2*temp2];
        k1 = k1+1;
    end
end
for i = 1:m-1
    for j = 1:n-1
        temp1 = nodes{i,j}(2,:) - nodes{i,j}(3,:);
        temp2 = nodes{i,j+1}(2,:) - nodes{i,j}(3,:);
        alpha1 = d/norm(temp1);
        alpha2 = d/norm(temp2);
        C1(k1,:) = [nodes{i,j}(2,:)-2*alpha1*temp1,nodes{i,j}(3,:),nodes{i,j+1}(2,:)-2*alpha2*temp2];
        k1 = k1+1;
    end
end
for i = 1:m-1
    for j = 1:n
        temp1 = nodes{i,j}(1,:) - nodes{i,j}(2,:);
        temp2 = nodes{i+1,j}(1,:) - nodes{i,j}(2,:);
        alpha1 = d/norm(temp1);
        alpha2 = d/norm(temp2);
        C1(k1,:) = [nodes{i,j}(1,:)-2*alpha1*temp1,nodes{i,j}(2,:),nodes{i+1,j}(1,:)-2*alpha2*temp2];
        k1 = k1+1;
    end
end
for i = 1:m
    for j = 1:n-1
        temp1 = nodes{i,j}(3,:) - nodes{i,j}(4,:);
        temp2 = nodes{i,j}(6,:) - nodes{i,j}(4,:);
        alpha1 = d/norm(temp1);
        alpha2 = d/norm(temp2);
        C1(k1,:) = [nodes{i,j}(3,:)-2*alpha1*temp1,nodes{i,j}(4,:),nodes{i,j}(6,:)-2*alpha2*temp2];
        k1 = k1+1;
    end
end
%build the half cut
for j = 1:n
    temp = nodes{1,j}(1,:) - nodes{1,j}(5,:);
    alpha = d/norm(temp);
    C2(k2,:) = [nodes{1,j}(1,:)-2*alpha*temp,nodes{1,j}(5,:)-2*alpha*temp];
    k2 = k2+1;
end
for j = 1:n
    temp = nodes{m,j}(2,:) - nodes{m,j}(1,:);
    alpha = d/norm(temp);
    C2(k2,:) = [nodes{m,j}(2,:)+2*alpha*temp,nodes{m,j}(1,:)+2*alpha*temp];
    k2 = k2+1;
end
for i = 1:m-1
    temp = nodes{i,1}(9,:) - nodes{i,1}(2,:);
    alpha = d/norm(temp);
    C2(k2,:) = [nodes{i,1}(9,:)+2*alpha*temp,nodes{i,1}(2,:)+2*alpha*temp];
    k2 = k2+1;
end
for i = 1:m-1
    temp = nodes{i,n}(3,:) - nodes{i,n}(2,:);
    alpha = d/norm(temp);
    C2(k2,:) = [nodes{i,n}(3,:)+2*alpha*temp,nodes{i,n}(2,:)+2*alpha*temp];
    k2 = k2+1;
end
%% build the geometry model in comsol
import com.comsol.model.*
import com.comsol.model.util.*
model = ModelUtil.create('Model');
model.modelPath(pwd);
model.component.create('comp1', true);
model.component('comp1').geom.create('geom1', 2);
model.component('comp1').mesh.create('mesh1');
model.component('comp1').physics.create('solid', 'SolidMechanics', 'geom1');
% set globle parameters
dis = num2str(displacement);
model.param.set('d', dis, 'the displacement in x direction');
model.component('comp1').geom('geom1').run;
% build the kirigami sheet
model.component('comp1').geom('geom1').create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('r1').set('size', [length_x length_y]);
model.component('comp1').geom('geom1').feature('r1').set('base', 'center');
model.component('comp1').geom('geom1').run('r1');
% build the cuts
for i = 1:k1-1
    fil = strcat('fil',num2str(i));
    thi = strcat('thi',num2str(i));
    ic = strcat('ic',num2str(i));
    model.component('comp1').geom('geom1').create(ic, 'InterpolationCurve');
    model.component('comp1').geom('geom1').feature(ic).set('table', [[linspace(C1(i,1),C1(i,3),10)';linspace(C1(i,3),C1(i,5),10)'],[linspace(C1(i,2),C1(i,4),10)';linspace(C1(i,4),C1(i,6),10)']]);
%     model.component('comp1').geom('geom1').feature(ic).set('table', [[inter_node(C1(i,1),C1(i,3),3,1/5)';flip(inter_node(C1(i,5),C1(i,3),3,1/5))'],[inter_node(C1(i,2),C1(i,4),3,1/5)';flip(inter_node(C1(i,6),C1(i,4),3,1/5))']]);
    model.component('comp1').geom('geom1').create(thi, 'Thicken2D');
    model.component('comp1').geom('geom1').feature(thi).set('totalthick', w);
    model.component('comp1').geom('geom1').feature(thi).set('ends', 'circular');
    model.component('comp1').geom('geom1').feature(thi).selection('input').set({ic});
%     model.component('comp1').geom('geom1').create(fil, 'Fillet');
%     model.component('comp1').geom('geom1').feature(fil).set('radius', w/2);
%     model.component('comp1').geom('geom1').feature(fil).selection('point').set(thi, [1 2 3 4]);
end
for i = 1:k2-1
    ls = strcat('ls',num2str(i));
    thi = strcat('thi',num2str(i+k1-1));
    fil = strcat('fil',num2str(i+k1-1));
    model.component('comp1').geom('geom1').create(ls, 'LineSegment');
    model.component('comp1').geom('geom1').feature(ls).set('specify1', 'coord');
    model.component('comp1').geom('geom1').feature(ls).set('coord1', C2(i,1:2));
    model.component('comp1').geom('geom1').feature(ls).set('specify2', 'coord');
    model.component('comp1').geom('geom1').feature(ls).set('coord2', C2(i,3:4));
    model.component('comp1').geom('geom1').create(thi, 'Thicken2D');
    model.component('comp1').geom('geom1').feature(thi).set('totalthick', w);
    model.component('comp1').geom('geom1').feature(thi).set('ends', 'circular');
    model.component('comp1').geom('geom1').feature(thi).selection('input').set({ls});
%     model.component('comp1').geom('geom1').create(fil, 'Fillet');
%     model.component('comp1').geom('geom1').feature(fil).set('radius', w/2);
%     model.component('comp1').geom('geom1').feature(fil).selection('point').set(thi, [1 2 3 4]);
end
model.component('comp1').geom('geom1').run;
% remove the cuts from kirigami sheet
model.component('comp1').geom('geom1').create('dif1', 'Difference');
model.component('comp1').geom('geom1').feature('dif1').selection('input').set({'r1'});
for i = 1:k1+k2-2
    thi = strcat('thi',num2str(i));
    model.component('comp1').geom('geom1').feature('dif1').selection('input2').add({thi});
end
model.component('comp1').geom('geom1').run('dif1');
model.component('comp1').geom('geom1').run;

% set material properties
model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material('mat1').propertyGroup.create('Enu', 'Young''s modulus and Poisson''s ratio');
model.component('comp1').material('mat1').propertyGroup('def').set('density', num2str(properties(1)));
model.component('comp1').material('mat1').propertyGroup('Enu').set('E', num2str(properties(2)));
model.component('comp1').material('mat1').propertyGroup('Enu').set('nu', num2str(properties(3)));
model.component('comp1').geom('geom1').run('fin');
% plot the geometry pattern
%figure(1)
%mphgeom(model) % plot the geometry
%% set the boundary conditions
model.component('comp1').physics('solid').create('rms1', 'RigidMotionSuppression', 2);
% set the displacement of the nodes on the right side
e0 = 1e-8;
coordBox1 = [n/2-e0 n/2+e0;-m/2-e0 m/2+e0];
index1 = mphselectbox(model,'geom1',coordBox1,'point');
coordBox2 = [-n/2-e0 -n/2+e0;-m/2-e0 m/2+e0];
index2 = mphselectbox(model,'geom1',coordBox2,'point');
[~,temp] = size(index);
% set the displacement of the nodes on the right side
model.component('comp1').physics('solid').create('disp1', 'Displacement0', 0);
model.component('comp1').physics('solid').feature('disp1').selection.set(index1);
model.component('comp1').physics('solid').feature('disp1').setIndex('Direction', true, 0);
model.component('comp1').physics('solid').feature('disp1').setIndex('U0', 'd', 0);
% set the displacement of the nodes on the left side
model.component('comp1').physics('solid').create('disp2', 'Displacement0', 0);
model.component('comp1').physics('solid').feature('disp2').selection.set(index2);
model.component('comp1').physics('solid').feature('disp2').setIndex('Direction', true, 0);
model.component('comp1').physics('solid').feature('disp2').setIndex('U0', '-d', 0);

% apply the thickness
model.component('comp1').physics('solid').prop('d').set('d', 0.003);
% mesh
model.component('comp1').mesh('mesh1').autoMeshSize(6);
model.component('comp1').mesh('mesh1').run;
% set the study
model.study.create('std1');
model.study('std1').create('stat', 'Stationary');
model.study('std1').feature('stat').set('geometricNonlinearity', true);

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature.remove('fcDef');

model.result.create('pg1', 'PlotGroup2D');
model.result('pg1').create('surf1', 'Surface');
model.result('pg1').feature('surf1').set('expr', 'solid.mises');
model.result('pg1').feature('surf1').create('def', 'Deform');

model.sol('sol1').attach('std1');
model.sol('sol1').feature('st1').label('Compile Equations: Stationary');
model.sol('sol1').feature('v1').label('Dependent Variables 1.1');
model.sol('sol1').feature('s1').label('Stationary Solver 1.1');
model.sol('sol1').feature('s1').feature('dDef').label('Direct 1');
model.sol('sol1').feature('s1').feature('aDef').label('Advanced 1');
model.sol('sol1').feature('s1').feature('aDef').set('cachepattern', true);
model.sol('sol1').feature('s1').feature('fc1').label('Fully Coupled 1.1');
model.sol('sol1').runAll;

model.result('pg1').label('Stress (solid)');
model.result('pg1').set('frametype', 'spatial');
model.result('pg1').feature('surf1').set('const', {'solid.refpntx' '0' 'Reference point for moment computation, x-coordinate'; 'solid.refpnty' '0' 'Reference point for moment computation, y-coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z-coordinate'});
model.result('pg1').feature('surf1').set('colortable', 'Prism');
model.result('pg1').feature('surf1').set('threshold', 'manual');
model.result('pg1').feature('surf1').set('thresholdvalue', 0.2);
model.result('pg1').feature('surf1').set('resolution', 'normal');
model.result('pg1').feature('surf1').feature('def').set('scaleactive', true);
model.sol('sol1').runAll;
%
%% plot the results
%plot the deformation graph
figure(2)
mphplot(model,'pg1','rangenum',1)
%select the bounday nodes
e0 = 1e-8;
coordBox = [-n/2-e0 n/2+e0;-m/2-e0 -m/2+e0];
index = mphselectbox(model,'geom1',coordBox,'point');
% save the boundary nodes data
x = mphevalpoint(model,'x','selection',index);   %2x4
y = mphevalpoint(model,'y','selection',index);
boundary_shape = [x' y'];
[a,~] = size(boundary_shape);
res = 0;
r = 3.2;
% distance control
% for i = 1:a
%     res1 = res + abs(1/sqrt(3)*sqrt(r^2-boundary_shape(i,1)^2) -  boundary_shape(i,2));
% end
% differential control
% for i = 1:a-1
%     vec1 = -1/sqrt(4)*sqrt(r^2-boundary_shape(i+1,1)^2) + 1/sqrt(4)*sqrt(r^2-boundary_shape(1,1)^2); % the distance of target shape
%     vec2 = boundary_shape(i+1,2) - boundary_shape(1,2);                                             % the distance of optimised shape
%     res = abs(vec2 - vec1) + res;
% end
% fit to a curve and calculate the difference area
x0 = boundary_shape(:,1);
y0 = boundary_shape(:,2);
p = polyfit(x0,y0,9); % fit the a ploynomial function
q = polyint(p);
% fun2 = 0;
% [~,iter] = size(q);
% for i=1:iter
%     fun2 = q(iter-i+1)*(sym(x).^(i-1))+fun2;
% end
%fun1 = @(x1) -1/sqrt(4)*sqrt(r^2-x1.^2);
%fun2 = @(x1) poly2sym(q,x1);
syms x1
fun = abs(-1/sqrt(6)*sqrt(r^2-x1^2)-poly2sym(q,x1));
res = int(fun,-m,m);

% curvature control
% fit the boundary nodes
boundary_shape = [x y];
