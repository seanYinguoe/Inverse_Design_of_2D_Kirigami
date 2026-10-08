function [res boundary_shape] = livelink_rigid(m,n,tessellation,displacement,properties,s)
% livelink matlab commands for comsol
%
% This model demonstrates how to use matlab optimize the kirigami pattern
% using the Livelink for matlab.
% units: m,
%% build geometry model in tessellation_optimization.m
temp = cell(m/2,n/2);

% derive the optimised variables
index = [];

% x coordinate
for j = 1:(n/2)
    index = horzcat(index,[5+(j-1)*10]);
    index = horzcat(index,[2+(j-1)*10+(n/2)*10*(m/2-1)]);
end
if n > 2
    for j = 1:(n/2-1)
        index = horzcat(index,[7+(j-1)*10]);
        index = horzcat(index,[3+(j-1)*10+(n/2)*10*(m/2-1)]);
    end
end

% y coordinate
for i = 1:(m/2)
    index = horzcat(index,[9+(i-1)*10*n/2+m*n/4*10]);
    index = horzcat(index,[4+(i-1)*10*n/2+10*(n/2-1)+m*n/4*10]);
end
if m > 2
    for i = 1:(m/2-1)
        index = horzcat(index,[10+(i-1)*10+m*n/4*10]);
        index = horzcat(index,[3+(i-1)*10+10*(n/2-1)+m*n/4*10]);
    end
end
% nodes extension from num_v variables to m*n/4*10*2 variables
x = zeros(m*n*20/4,1);
x(index) = tessellation(1:8);
for i = 1:(m/2)
    x([8,9,10]+(n/2)*(i-1)*10) = -n/2;
end
for i = 1:m/2
    for j = 1:(n/2-1)
        x([18 20]+10*(j-1)+(i-1)*(n/2)*10) = x([7 3]+10*(j-1)+(i-1)*(n/2)*10);
        x([18 20]+10*(j-1)+10*m*n/4+(i-1)*(n/2)*10) = x([7 3]+10*(j-1)+10*m*n/4+(i-1)*(n/2)*10);
    end
end
for j = 1:n/2
    x([5 7 8]+10*(j-1)+10*m*n/4) = -m/2;
end
% add nodes 4, 9, 1, 6
node_4 = tessellation(9) * ([tessellation(6),0]-[tessellation(5),-m/2]) + [tessellation(5),-m/2];
node_9 = tessellation(10) * ([tessellation(6),0]-[tessellation(5),-m/2]) + [tessellation(5),-m/2];
node_1 = [tessellation(2) (tessellation(2)+n/2)/(node_4(1)+n/2)*(node_4(2)-tessellation(7))+tessellation(7)];
node_6 = tessellation(11) * (node_4 - [-n/2 tessellation(7)]) + [-n/2 tessellation(7)];
node_21 = [tessellation(4) (tessellation(4)-node_9(1))/(0-node_9(1))*(tessellation(8)-node_9(2))+node_9(2)];
node_26 = tessellation(12) * ([0 tessellation(8)] - node_9) + node_9;
x([1 21 4 24 19 39 6 26 16 36 11 31]) = [node_1 node_4 node_9 node_6 node_26 node_21];

% derive the nodes (m/2)*(n/2)
for i = 1:m/2
    for j = 1:n/2
        temp{i,j}(:,1) = [x((i-1)*n/2*10+10*(j-1)+1:(i-1)*n/2*10+10*j)];
        temp{i,j}(:,2) = [x(10*m*n/4+(i-1)*n/2*10+10*(j-1)+1:10*m*n/4+(i-1)*n/2*10+10*j)];
    end
end
nodes = derive_nodes_rigid(temp); % dupicate nodes (m/2)*(n/2) to m*n
length_x = nodes{1,n}(7,1) - nodes{1,1}(8,1);
length_y = nodes{m,1}(10,2) - nodes{1,1}(8,2);
d = 0.035; %the gap between cuts 0.03
w = 0.02; %the width of cuts
f_nodes = [1 2 3 4];                % set the fixed nodes
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
%% build the geometry model in comsol
import com.comsol.model.*
import com.comsol.model.util.*
model = ModelUtil.create('Model');
model.modelPath(pwd);
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
    ls = strcat('ls',num2str(i));
    thi = strcat('thi',num2str(i));
    fil = strcat('fil',num2str(i));
    model.component('comp1').geom('geom1').create(ls, 'LineSegment');
    model.component('comp1').geom('geom1').feature(ls).set('specify1', 'coord');
    model.component('comp1').geom('geom1').feature(ls).set('coord1', C1(i,1:2));
    model.component('comp1').geom('geom1').feature(ls).set('specify2', 'coord');
    model.component('comp1').geom('geom1').feature(ls).set('coord2', C1(i,3:4));
    model.component('comp1').geom('geom1').create(thi, 'Thicken2D');
    model.component('comp1').geom('geom1').feature(thi).set('totalthick', w);
    model.component('comp1').geom('geom1').feature(thi).selection('input').set({ls});
    model.component('comp1').geom('geom1').create(fil, 'Fillet');
    model.component('comp1').geom('geom1').feature(fil).set('radius', w/4);
    model.component('comp1').geom('geom1').feature(fil).selection('point').set(thi, [1 2 3 4]);
end
% build the half cuts
for i = 1:k2-1
    ls = strcat('ls',num2str(i+k1-1));
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
    model.component('comp1').geom('geom1').create(fil, 'Fillet');
    model.component('comp1').geom('geom1').feature(fil).set('radius', w/4);
    model.component('comp1').geom('geom1').feature(fil).selection('point').set(thi, [1 2 3 4]);
end
model.component('comp1').geom('geom1').run;
% remove the cuts from kirigami sheet
model.component('comp1').geom('geom1').create('dif1', 'Difference');
model.component('comp1').geom('geom1').feature('dif1').selection('input').set({'r1'});
for i = 1:k1+k2-2
%     thi = strcat('thi',num2str(i));
%     model.component('comp1').geom('geom1').feature('dif1').selection('input2').add({thi});
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
if s== 1
    % convex
    % Half cut bending energy
    energy_total = 0;
    select_nodes1 = mphgetcoords(model,'geom1','point',[34 30 40 32 94 90 100 92 14 3 13 2 33 29 39 31 93 89 99 91]); % extract specified node coordinates from comsol 1/2
    [~,count1] = size(select_nodes1);
    for i = 1:count1
        ux = mphinterp(model,'u','coord',select_nodes1(:,i),'dataset', 'dset1'); % extract displacement in x direction
        uy = mphinterp(model,'v','coord',select_nodes1(:,i),'dataset', 'dset1'); % extract displacement in y direction
        nodes1{i} = select_nodes1(:,i) + [ux'; uy'];
    end
    for i = 1:count1/4
        theta1{i} = bending_angle(nodes1(4*(i-1)+1:4*i)); % get bending angle of each hinge
    end
    for i = 1:size(theta1,2)
        energy = 2 * bending_energy(theta1{i}-theta1{i}(1));  % get bending energy of each hinge
        energy_total = energy_total + energy;
    end
    % Get the whole cut bending energy(substract the initial angle)
    select_nodes2 = mphgetcoords(model,'geom1','point',[10 36 12 20 58 38 60 22 66 62 72 50 42 64 48 52 54 96 56 80 102 98 104 82]); % extract 1/4 nodes
    select_nodes3 = mphgetcoords(model,'geom1','point',[18 16 24 26 17 15 23 25 28 44 27 43 74 46 73 45 78 76 84 86 77 75 83 85 112 106 124 128 111 108 123 126]); % extract 1/2 nodes
    select_nodes4 = mphgetcoords(model,'geom1','point',[88 115 87 114 146 119 145 118]); % extract nodes
    [~,count2] = size(select_nodes2);
    [~,count3] = size(select_nodes3);
    [~,count4] = size(select_nodes4);
    % Get the 1/4 whole cut energy
    for i = 1:count2
        ux = mphinterp(model,'u','coord',select_nodes2(:,i),'dataset', 'dset1'); % extract displacement in x direction
        uy = mphinterp(model,'v','coord',select_nodes2(:,i),'dataset', 'dset1'); % extract displacement in y direction
        nodes2{i} = select_nodes2(:,i) + [ux'; uy'];
    end
    for i = 1:count2/4
        theta2{i} = bending_angle(nodes2(4*(i-1)+1:4*i)); % get bending angle of each hinge
    end
    for i = 1:size(theta2,2)
        energy = 4 * bending_energy(theta2{i}-theta2{i}(1));  % get bending energy of each hinge
        energy_total = energy_total + energy;
    end
    % Get the 1/2 whole cut energy
    for i = 1:count3
        ux = mphinterp(model,'u','coord',select_nodes3(:,i),'dataset', 'dset1'); % extract displacement in x direction
        uy = mphinterp(model,'v','coord',select_nodes3(:,i),'dataset', 'dset1'); % extract displacement in y direction
        nodes3{i} = select_nodes3(:,i) + [ux'; uy'];
    end
    for i = 1:count3/4
        theta3{i} = bending_angle(nodes3(4*(i-1)+1:4*i)); % get bending angle of each hinge
    end
    for i = 1:size(theta3,2)
        energy = 2 * bending_energy(theta3{i}-theta3{i}(1));  % get bending energy of each hinge
        energy_total = energy_total + energy;
    end
    % Get the whole cut energy
    for i = 1:count4
        ux = mphinterp(model,'u','coord',select_nodes4(:,i),'dataset', 'dset1'); % extract displacement in x direction
        uy = mphinterp(model,'v','coord',select_nodes4(:,i),'dataset', 'dset1'); % extract displacement in y direction
        nodes4{i} = select_nodes4(:,i) + [ux'; uy'];
    end
    for i = 1:count4/4
        theta4{i} = bending_angle(nodes4(4*(i-1)+1:4*i)); % get bending angle of each hinge
    end
    for i = 1:size(theta4,2)
        energy = bending_energy(theta4{i}-theta4{i}(1));  % get bending energy of each hinge
        energy_total = energy_total + energy;
    end
elseif s==2
    % concave
    % Half cut bending energy
    energy_total = 0;
    select_nodes1 = mphgetcoords(model,'geom1','point',[30 46 36 48 90 98 96 100 14 3 13 2 45 29 47 35 97 89 99 95]); % extract specified node coordinates from comsol 1/2
    [~,count1] = size(select_nodes1);
    for i = 1:count1
        ux = mphinterp(model,'u','coord',select_nodes1(:,i),'dataset', 'dset1'); % extract displacement in x direction
        uy = mphinterp(model,'v','coord',select_nodes1(:,i),'dataset', 'dset1'); % extract displacement in y direction
        nodes1{i} = select_nodes1(:,i) + [ux'; uy'];
    end
    for i = 1:count1/4
        theta1{i} = bending_angle(nodes1(4*(i-1)+1:4*i)); % get bending angle of each hinge
    end
    for i = 1:size(theta1,2)
        energy = 2 * bending_energy(theta1{i}-theta1{i}(1));  % get bending energy of each hinge
        energy_total = energy_total + energy;
    end
    % Get the whole cut bending energy(substract the initial angle)
    select_nodes2 = mphgetcoords(model,'geom1','point',[10 32 12 20 58 34 60 22 82 62 88 50 38 64 44 52 54 92 56 72 102 94 104 74]); % extract 1/4 nodes
    select_nodes3 = mphgetcoords(model,'geom1','point',[18 16 24 26 17 15 23 25 28 40 27 39 66 42 65 41 70 68 76 78 69 67 75 77 112 106 124 128 111 108 123 126]); % extract 1/2 nodes
    select_nodes4 = mphgetcoords(model,'geom1','point',[80 115 79 114 154 119 153 118]); % extract nodes
    [~,count2] = size(select_nodes2);
    [~,count3] = size(select_nodes3);
    [~,count4] = size(select_nodes4);
    % Get the 1/4 whole cut energy
    for i = 1:count2
        ux = mphinterp(model,'u','coord',select_nodes2(:,i),'dataset', 'dset1'); % extract displacement in x direction
        uy = mphinterp(model,'v','coord',select_nodes2(:,i),'dataset', 'dset1'); % extract displacement in y direction
        nodes2{i} = select_nodes2(:,i) + [ux'; uy'];
    end
    for i = 1:count2/4
        theta2{i} = bending_angle(nodes2(4*(i-1)+1:4*i)); % get bending angle of each hinge
    end
    for i = 1:size(theta2,2)
        energy = 4 * bending_energy(theta2{i}-theta2{i}(1));  % get bending energy of each hinge
        energy_total = energy_total + energy;
    end
    % Get the 1/2 whole cut energy
    for i = 1:count3
        ux = mphinterp(model,'u','coord',select_nodes3(:,i),'dataset', 'dset1'); % extract displacement in x direction
        uy = mphinterp(model,'v','coord',select_nodes3(:,i),'dataset', 'dset1'); % extract displacement in y direction
        nodes3{i} = select_nodes3(:,i) + [ux'; uy'];
    end
    for i = 1:count3/4
        theta3{i} = bending_angle(nodes3(4*(i-1)+1:4*i)); % get bending angle of each hinge
    end
    for i = 1:size(theta3,2)
        energy = 2 * bending_energy(theta3{i}-theta3{i}(1));  % get bending energy of each hinge
        energy_total = energy_total + energy;
    end
    % Get the whole cut energy
    for i = 1:count4
        ux = mphinterp(model,'u','coord',select_nodes4(:,i),'dataset', 'dset1'); % extract displacement in x direction
        uy = mphinterp(model,'v','coord',select_nodes4(:,i),'dataset', 'dset1'); % extract displacement in y direction
        nodes4{i} = select_nodes4(:,i) + [ux'; uy'];
    end
    for i = 1:count4/4
        theta4{i} = bending_angle(nodes4(4*(i-1)+1:4*i)); % get bending angle of each hinge
    end
    for i = 1:size(theta4,2)
        energy = bending_energy(theta4{i}-theta4{i}(1));  % get bending energy of each hinge
        energy_total = energy_total + energy;
    end
end
%% plot the results
%plot the deformation graph
%figure(2);
%mphplot(model,'pg1','rangenum',1);
%select the bounday nodes
e0 = 1e-8;
if s == 3
    coordBox = [-n/2-e0 n/2+e0;-m/2-e0 -m/2+e0];
elseif s == 1 | s == 2
    coordBox = [-n/2+e0 n/2-e0;-m/2-e0 -m/2+e0];
end
index = mphselectbox(model,'geom1',coordBox,'point');
% save the boundary nodes data
x = mphevalpoint(model,'x','selection',index);
y = mphevalpoint(model,'y','selection',index);
x = x(:,end);
y = y(:,end);
boundary_shape = [x' y'];
[~,a] = size(boundary_shape);
res = 0;
r = 3.2*sqrt(1);

%% calculate the distance between target shape and boudary shape

%% fit to a curve and calculate the difference area
x0 = x;
y0 = y;
p = polyfit(x0,y0,5); % fit the a ploynomial function: ellipse 7; sinwave; vase 5
%% calculate the error
res = shape_difference(p,s,r)
%res = shape_difference(s,r,boundary_shape)
end
