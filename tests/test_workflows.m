function tests = test_workflows
tests=functiontests(localfunctions);
end
function setupOnce(t)
t.applyFixture(matlab.unittest.fixtures.PathFixture(fileparts(fileparts(mfilename('fullpath')))));setup_project;
end
function testPresetAndNumericBoundaryAgree(t)
T=tessellation_deployment(2,4,1,.5);
for code=1:5
    target=make_target(code,'Size',3.2);
    for mode={'rigid','nonrigid'}
        [c,e]=boundary_residual(T,code,3.2,2,4,mode{1});
        [cc,ee]=boundary_residual(T,target,3.2,2,4,mode{1});
        t.verifyEqual(cc,c);t.verifyEqual(ee,e);
    end
end
end
function testWavyEquationAndProfileAgree(t)
a=make_target('wavy');x=linspace(-2.5,2.5,51)';y=target_profile(a,x,'top');
t.verifyEqual(shape(4,[x y],a.r),zeros(size(x)),'AbsTol',1e-12);
end
function testMechanicalRoundTrip(t)
for mode={'rigid','nonrigid'}
    k=seed(2,4,mode{1});p=prepare_mechanical_pattern(k);
    t.verifyLessThan(p.reconstruction_error,1e-12);
    nodes=decode_pattern(p.x0,p);
    if strcmp(mode{1},'rigid'),[a,b]=cuts_rigid(nodes,.035);else,[a,b]=cuts_nonrigid(nodes,.03);end
    t.verifySize(a,[21 4+2*strcmp(mode{1},'nonrigid')]);t.verifySize(b,[10 4]);
end
end
function testNonrigidRectangularIndexing(t)
k=seed(4,6,'nonrigid');k.config.target=make_target('vase','Size',4.2,'ClampX',3.5);
p=prepare_mechanical_pattern(k);t.verifyLessThan(p.reconstruction_error,1e-12);
t.verifyEqual(numel(p.x0),2*4*6-2);
end
function testRejectUnrepresentableSeed(t)
k=seed(2,4,'nonrigid');k.compact{2,4}(1,1)=k.compact{2,4}(1,1)+.1;
t.verifyError(@()prepare_mechanical_pattern(k),'kirigami:UnrepresentableSeed');
end
function testShapeErrorIgnoresVerticalTranslation(t)
target=make_target('ellipse','Size',3.2);x=linspace(-2.5,2.5,21)';y=.1*x.^2-1.5;
[e,~]=profile_error([x y],target,2);[f,~]=profile_error([x y+4],target,2);
t.verifyEqual(e,f,'AbsTol',1e-10);
t.verifyGreaterThan(e,0);
end
function testGAElitismAndBounds(t)
rng(3);cfg=mechanical_config;g=cfg.ga;g.population=8;g.generations=4;g.tolerance=0;
r=genetic_search(@quadratic,[.8 -.6],[-1 -1],[1 1],g);
t.verifyLessThanOrEqual(diff(r.history(:,2)),zeros(size(r.history,1)-1,1));
t.verifyTrue(all(r.x>=-1&r.x<=1));t.verifyLessThanOrEqual(r.fitness,1);
end
function testGAZeroFitness(t)
c=mechanical_config;c.ga.population=3;
r=genetic_search(@quadratic,[0 0],[-1 -1],[1 1],c.ga);
t.verifyEqual(r.fitness,0);t.verifyEqual(r.generation,0);t.verifyEqual(r.stop_reason,'target tolerance');
end
function testGAFailedBaselineIsNotHidden(t)
c=mechanical_config;c.ga.population=3;
t.verifyError(@()genetic_search(@fail,[0],[-1],[1],c.ga),'test:Baseline');
end
function [f,d]=quadratic(x)
f=sum(x.^2);d=struct('x',x);
end
function [f,d]=fail(~) %#ok<STOUT>
error('test:Baseline','Baseline must fail visibly.');
end
function k=seed(m,n,mode)
cfg=kinematic_config;cfg.mode=mode;
k=struct('compact',{tessellation_deployment(m,n,1,0)}, ...
    'config',cfg,'diagnostics',struct('converged',true));
end
