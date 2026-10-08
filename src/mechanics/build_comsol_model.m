function model = build_comsol_model(C1,C2,pattern,cfg)
%BUILD_COMSOL_MODEL Build the source 2-D nonlinear solid model entirely in MATLAB.
% No .mph input is needed and no .mph output is written.
import com.comsol.model.util.*
validateattributes(cfg.metres_per_unit,{'numeric'},{'scalar','finite','positive'});
validateattributes(cfg.load_steps,{'numeric'},{'scalar','integer','positive'});
validateattributes(cfg.thickness,{'numeric'},{'scalar','finite','positive'});
validateattributes(cfg.cut_width,{'numeric'},{'scalar','finite','positive'});
validateattributes(cfg.material.young_modulus,{'numeric'},{'scalar','finite','positive'});
validateattributes(cfg.material.poisson_ratio,{'numeric'},{'scalar','>',-1,'<',0.5});
scale=cfg.metres_per_unit;C1=C1*scale;C2=C2*scale;w=cfg.cut_width*scale;
m=pattern.rows;n=pattern.cols;
properties=[cfg.material.density cfg.material.young_modulus cfg.material.poisson_ratio];
ramp=strtrim(sprintf('%.16g ',linspace(0,cfg.displacement*scale,cfg.load_steps+1)));
tag=char(ModelUtil.uniquetag('Kirigami'));model=ModelUtil.create(tag);
try
model.hist.disable;
model.component.create('comp1',true);
model.component('comp1').geom.create('geom1',2);
model.component('comp1').geom('geom1').lengthUnit('m');
model.component('comp1').mesh.create('mesh1');
model.component('comp1').physics.create('solid','SolidMechanics','geom1');
model.param.set('d','0','Per-grip displacement in metres');
g=model.component('comp1').geom('geom1');
g.create('r1','Rectangle');g.feature('r1').set('size',[n m]*scale);
g.feature('r1').set('base','center');
paths=[mat2cell(C1,ones(size(C1,1),1),size(C1,2)); ...
    mat2cell(C2,ones(size(C2,1),1),4)];
for k=1:numel(paths)
    p=reshape(paths{k},2,[])';
    assert(all(vecnorm(diff(p),2,2)>w),'A cut segment is shorter than its width.');
    line=sprintf('cut%d',k);thick=sprintf('thick%d',k);fil=sprintf('fil%d',k);
    if size(p,1)==2
        g.create(line,'LineSegment');
        g.feature(line).set('specify1','coord');g.feature(line).set('coord1',p(1,:));
        g.feature(line).set('specify2','coord');g.feature(line).set('coord2',p(2,:));
    else
        g.create(line,'InterpolationCurve');
        q=linspace(0,1,10)';
        table=[p(1,:)+q.*(p(2,:)-p(1,:));p(2,:)+q(2:end).*(p(3,:)-p(2,:))];
        g.feature(line).set('table',table);
    end
    g.create(thick,'Thicken2D');g.feature(thick).set('totalthick',w);
    g.feature(thick).selection('input').set({line});
    g.create(fil,'Fillet');
    if strcmp(pattern.mode,'rigid'),radius=w/4;else,radius=w/8;end
    g.feature(fil).set('radius',radius);
    g.feature(fil).selection('point').set(thick,[1 2 3 4]);
end
g.create('dif1','Difference');g.feature('dif1').selection('input').set({'r1'});
for k=1:numel(paths),g.feature('dif1').selection('input2').add({sprintf('fil%d',k)});end
g.run;
assert(g.getNDomains()==1,'kirigami:DisconnectedSheet', ...
    'Cut geometry must leave one connected sheet.');
model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material('mat1').propertyGroup.create('Enu', 'Young''s modulus and Poisson''s ratio');
model.component('comp1').material('mat1').propertyGroup('def').set('density', num2str(properties(1)));
model.component('comp1').material('mat1').propertyGroup('Enu').set('E', num2str(properties(2)));
model.component('comp1').material('mat1').propertyGroup('Enu').set('nu', num2str(properties(3)));
model.component('comp1').geom('geom1').run('fin');
model.component('comp1').physics('solid').create('rms1', 'RigidMotionSuppression', 2);
e0 = 1e-8*scale;
coordBox1 = [n*scale/2-e0 n*scale/2+e0;-m*scale/2-e0 m*scale/2+e0];
index1 = mphselectbox(model,'geom1',coordBox1,'point');
coordBox2 = [-n*scale/2-e0 -n*scale/2+e0;-m*scale/2-e0 m*scale/2+e0];
index2 = mphselectbox(model,'geom1',coordBox2,'point');
assert(~isempty(index1)&&~isempty(index2),'kirigami:MissingGrips', ...
    'Could not find both loading grips.');

model.component('comp1').physics('solid').create('disp1', 'Displacement0', 0);
model.component('comp1').physics('solid').feature('disp1').selection.set(index1);
model.component('comp1').physics('solid').feature('disp1').setIndex('Direction', true, 0);
model.component('comp1').physics('solid').feature('disp1').setIndex('U0', 'd', 0);
model.component('comp1').physics('solid').create('disp2', 'Displacement0', 0);
model.component('comp1').physics('solid').feature('disp2').selection.set(index2);
model.component('comp1').physics('solid').feature('disp2').setIndex('Direction', true, 0);
model.component('comp1').physics('solid').feature('disp2').setIndex('U0', '-d', 0);
model.component('comp1').physics('solid').create('hmm1', 'HyperelasticModel', 2);
model.component('comp1').physics('solid').feature('hmm1').selection.all;

model.component('comp1').physics('solid').prop('d').set('d', cfg.thickness);
model.component('comp1').mesh('mesh1').autoMeshSize(cfg.mesh_size);
model.component('comp1').mesh('mesh1').run;
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


model.sol('sol1').attach('std1');
model.sol('sol1').feature('st1').label('Compile Equations: Stationary');
model.sol('sol1').feature('v1').label('Dependent Variables 1.1');
model.sol('sol1').feature('v1').set('clistctrl', {'p1'});
model.sol('sol1').feature('v1').set('cname', {'d'});
model.sol('sol1').feature('v1').set('clist', {ramp});
model.sol('sol1').feature('s1').label('Stationary Solver 1.1');
model.sol('sol1').feature('s1').set('probesel', 'none');
model.sol('sol1').feature('s1').feature('dDef').label('Direct 1');
model.sol('sol1').feature('s1').feature('aDef').label('Advanced 1');
model.sol('sol1').feature('s1').feature('aDef').set('cachepattern', true);
model.sol('sol1').feature('s1').feature('p1').label('Parametric 1.1');
model.sol('sol1').feature('s1').feature('p1').set('pname', {'d'});
model.sol('sol1').feature('s1').feature('p1').set('plistarr', {ramp});
model.sol('sol1').feature('s1').feature('p1').set('punit', {''});
model.sol('sol1').feature('s1').feature('p1').set('plot', false);
model.sol('sol1').feature('s1').feature('fc1').label('Fully Coupled 1.1');

catch e
    ModelUtil.remove(tag);rethrow(e);
end
end
