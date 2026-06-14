function [res boundary_shape] = livelink_new(m,n,tessellation,displacement,properties,s)
% livelink matlab commands for comsol
%
% This model demonstrates how to use matlab optimize the kirigami pattern 
% using the Livelink for matlab.
% units: m,
%% build geometry model in tessellation_optimization.m
temp = cell(m/2,n/2);
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
% nodes extension from num_v variables to m*n/4*9*2 variables
x = zeros(m*n*18/4,1);
x(index) = tessellation;
for i = 1:(m/2)
    x([7,8,9]+(n/2)*(i-1)*9) = -n/2;
end
%x([16,17,18,35]) = x([6,4,3,22]);
for i = 1:m/2
    for j = 1:(n/2-1)
        x([18 17 16]+9*(j-1)+(i-1)*(n/2)*9) = x([3 4 6]+9*(j-1)+(i-1)*(n/2)*9);
        x([18 17 16]+9*(j-1)+9*m*n/4+(i-1)*(n/2)*9) = x([3 4 6]+9*(j-1)+9*m*n/4+(i-1)*(n/2)*9);
    end
end
for i = 1:m/2-1
    for j = 1:n/2
        x([23 24 25]+9*(j-1)+(i-1)*(n/2)*9) = x([2 3 9]+9*(j-1)+(i-1)*(n/2)*9);
        x([23 24 25]+9*(j-1)+(i-1)*(n/2)*9+9*m*n/4) = x([2 3 9]+9*(j-1)+(i-1)*(n/2)*9+9*m*n/4);
    end
end
for j = 1:n/2
    x([5 6 7]+9*(j-1)+9*m*n/4) = -m/2;
end
%x([23,24,25,32,33,34]) = -1;
% derive the nodes (m/2)*(n/2)
for i = 1:m/2
    for j = 1:n/2
        temp{i,j}(:,1) = [x((i-1)*n/2*9+9*(j-1)+1:(i-1)*n/2*9+9*j)];
        temp{i,j}(:,2) = [x(9*m*n/4+(i-1)*n/2*9+9*(j-1)+1:9*m*n/4+(i-1)*n/2*9+9*j)]; 
    end
end
nodes = derive_nodes(temp); % dupicate nodes (m/2)*(n/2) to m*n
length_x = nodes{1,n}(6,1) - nodes{1,1}(7,1);
length_y = nodes{m,1}(9,2) - nodes{1,1}(7,2);
d = 0.02; %the gap between cuts 0.04
w = 0.02; %the width of cuts 0.02
f_nodes = [1 2 3 4];                % set the fixed nodes
%d_nodes = [177 178 179 180];        % set the displacement nodes 2*4
%d_nodes = [269 270 271 272];        %set the displacement nodes 2*6
% set the material properties
density = 2700;                 % density of material
E = 4.33e9;                       % the elastic module of material
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
        C1(k1,:) = [nodes{i,j}(8,:)-alpha1*temp1,nodes{i,j}(1,:),nodes{i,j}(4,:)-alpha2*temp2];
        k1 = k1+1;
    end
end
for i = 1:m-1
    for j = 1:n-1
        temp1 = nodes{i,j}(2,:) - nodes{i,j}(3,:);
        temp2 = nodes{i,j+1}(2,:) - nodes{i,j}(3,:);
        alpha1 = d/norm(temp1);
        alpha2 = d/norm(temp2);
        C1(k1,:) = [nodes{i,j}(2,:)-alpha1*temp1,nodes{i,j}(3,:),nodes{i,j+1}(2,:)-alpha2*temp2];
        k1 = k1+1;
    end
end
for i = 1:m-1
    for j = 1:n
        temp1 = nodes{i,j}(1,:) - nodes{i,j}(2,:);
        temp2 = nodes{i+1,j}(1,:) - nodes{i,j}(2,:);
        alpha1 = d/norm(temp1);
        alpha2 = d/norm(temp2);
        C1(k1,:) = [nodes{i,j}(1,:)-alpha1*temp1,nodes{i,j}(2,:),nodes{i+1,j}(1,:)-alpha2*temp2];
        k1 = k1+1;
    end
end
for i = 1:m
    for j = 1:n-1
        temp1 = nodes{i,j}(3,:) - nodes{i,j}(4,:);
        temp2 = nodes{i,j}(6,:) - nodes{i,j}(4,:);
        alpha1 = d/norm(temp1);
        alpha2 = d/norm(temp2);
        C1(k1,:) = [nodes{i,j}(3,:)-alpha1*temp1,nodes{i,j}(4,:),nodes{i,j}(6,:)-alpha2*temp2];
        k1 = k1+1;
    end
end
%build the half cut
for j = 1:n
    temp = nodes{1,j}(1,:) - nodes{1,j}(5,:);
    C2(k2,:) = [nodes{1,j}(1,:),nodes{1,j}(5,:)];
    k2 = k2+1;
end
for j = 1:n
    temp = nodes{m,j}(2,:) - nodes{m,j}(1,:);
    C2(k2,:) = [nodes{m,j}(2,:),nodes{m,j}(1,:)];
    k2 = k2+1;
end
for i = 1:m-1
    temp = nodes{i,1}(9,:) - nodes{i,1}(2,:);
    C2(k2,:) = [nodes{i,1}(9,:),nodes{i,1}(2,:)];
    k2 = k2+1;
end
for i = 1:m-1
    temp = nodes{i,n}(3,:) - nodes{i,n}(2,:);
    C2(k2,:) = [nodes{i,n}(3,:),nodes{i,n}(2,:)];
    k2 = k2+1;
end
%% build the geometry model in comsol
import com.comsol.model.*
import com.comsol.model.util.*
model = ModelUtil.create('Model');
model.modelPath('/Users/sean/Desktop/program/comsol');
model.hist.disable; % disable the model history
model.disableUpdates(true); % disable model feature update 
model.component.create('comp1', true);
model.component('comp1').geom.create('geom1', 2);
model.component('comp1').mesh.create('mesh1');
model.component('comp1').physics.create('solid', 'SolidMechanics', 'geom1');
% set globle parameters
dis = num2str(displacement);
model.param.set('d', '0', 'the displacement in x direction');
model.component('comp1').geom('geom1').run;
% build the kirigami sheet
model.component('comp1').geom('geom1').create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('r1').set('size', [length_x length_y]);
model.component('comp1').geom('geom1').feature('r1').set('base', 'center');
model.component('comp1').geom('geom1').run('r1');
% build the whole cuts
for i = 1:k1-1
    fil = strcat('fil',num2str(i));
    thi = strcat('thi',num2str(i));
    ic = strcat('ic',num2str(i));
    model.component('comp1').geom('geom1').create(ic, 'InterpolationCurve');
    model.component('comp1').geom('geom1').feature(ic).set('table', [[linspace(C1(i,1),C1(i,3),10)';linspace(C1(i,3),C1(i,5),10)'],[linspace(C1(i,2),C1(i,4),10)';linspace(C1(i,4),C1(i,6),10)']]);
    model.component('comp1').geom('geom1').create(thi, 'Thicken2D');
    model.component('comp1').geom('geom1').feature(thi).set('totalthick', w);
    model.component('comp1').geom('geom1').feature(thi).selection('input').set({ic});
    model.component('comp1').geom('geom1').create(fil, 'Fillet');
    model.component('comp1').geom('geom1').feature(fil).set('radius', w/8);
    model.component('comp1').geom('geom1').feature(fil).selection('point').set(thi, [1 2 3 4]);
end
% build the half cuts
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
    model.component('comp1').geom('geom1').feature(thi).selection('input').set({ls});
end
model.component('comp1').geom('geom1').run;
% remove the cuts from kirigami sheet
model.component('comp1').geom('geom1').create('dif1', 'Difference');
model.component('comp1').geom('geom1').feature('dif1').selection('input').set({'r1'});
for i = 1:k1+k2-2
    fil = strcat('fil',num2str(i));
    model.component('comp1').geom('geom1').feature('dif1').selection('input2').add({fil});
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
figure(1)
mphgeom(model) % plot the geometry
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
model.component('comp1').physics('solid').create('hmm1', 'HyperelasticModel', 2);
model.component('comp1').physics('solid').feature('hmm1').selection.set([1]);

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
model.sol('sol1').feature('s1').create('p1', 'Parametric');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature.remove('fcDef');

% model.result.create('pg1', 'PlotGroup2D');
% model.result('pg1').create('surf1', 'Surface');
% model.result('pg1').feature('surf1').set('expr', 'solid.mises');
% model.result('pg1').feature('surf1').create('def', 'Deform');

model.sol('sol1').attach('std1');
model.sol('sol1').feature('st1').label('Compile Equations: Stationary');
model.sol('sol1').feature('v1').label('Dependent Variables 1.1');
model.sol('sol1').feature('v1').set('clistctrl', {'p1'});
model.sol('sol1').feature('v1').set('cname', {'d'});
model.sol('sol1').feature('v1').set('clist', {'range(0,0.01,0.5)'});
model.sol('sol1').feature('s1').label('Stationary Solver 1.1');
model.sol('sol1').feature('s1').set('probesel', 'none');
model.sol('sol1').feature('s1').feature('dDef').label('Direct 1');
model.sol('sol1').feature('s1').feature('aDef').label('Advanced 1');
model.sol('sol1').feature('s1').feature('aDef').set('cachepattern', true);
model.sol('sol1').feature('s1').feature('p1').label('Parametric 1.1');
model.sol('sol1').feature('s1').feature('p1').set('pname', {'d'});
model.sol('sol1').feature('s1').feature('p1').set('plistarr', {'range(0,0.01,0.5)'});
model.sol('sol1').feature('s1').feature('p1').set('punit', {''});
model.sol('sol1').feature('s1').feature('p1').set('plot', true);
model.sol('sol1').feature('s1').feature('fc1').label('Fully Coupled 1.1');

model.result.create('pg1', 'PlotGroup2D');
model.result('pg1').set('data', 'dset1');
model.result('pg1').set('defaultPlotID', 'stress');
model.result('pg1').label('Stress (solid)');
model.result('pg1').set('frametype', 'spatial');
model.result('pg1').create('surf1', 'Surface');
model.result('pg1').feature('surf1').set('expr', {'solid.mises'});
model.result('pg1').feature('surf1').set('threshold', 'manual');
model.result('pg1').feature('surf1').set('thresholdvalue', 0.2);
model.result('pg1').feature('surf1').set('resolution', 'normal');
model.result('pg1').feature('surf1').set('colortable', 'Prism');
model.result('pg1').feature('surf1').create('def', 'Deform');
model.result('pg1').feature('surf1').feature('def').set('scaleactive', true);
model.result('pg1').feature('surf1').feature('def').set('scale', '1');
model.result('pg1').feature('surf1').feature('def').set('expr', {'u' 'v'});
model.result('pg1').feature('surf1').feature('def').set('descr', 'Displacement field');
model.sol('sol1').runAll;
% 
%% plot the results
%plot the deformation graph
figure(2);
mphplot(model,'pg1','rangenum',1);
%select the bounday nodes
e0 = 1e-8;
if s == 1 | s == 3
    coordBox = [-n/2-e0 n/2+e0;-m/2-e0 -m/2+e0];
elseif s == 2
    coordBox = [-n/2+e0 n/2-e0;-m/2-e0 -m/2+e0];
end
index = mphselectbox(model,'geom1',coordBox,'point');
% save the boundary nodes data
x = mphevalpoint(model,'x','selection',index);   
y = mphevalpoint(model,'y','selection',index);
x = x(:,end);
y = y(:,end);
boundary_shape = [x' y'];
[a,~] = size(boundary_shape);
res = 0;
r = 3.2*sqrt(0.9);
%% fit to a curve and calculate the difference area
x0 = x;
y0 = y;
p = polyfit(x0,y0,6); % fit the a ploynomial function: ellipse 7; sinwave 
%% calculate the error
res = shape_difference(p,s,r)
end

