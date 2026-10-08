function result = main_comsol(input,cfg)
%MAIN_COMSOL Refine a saved kinematic result, or the result struct itself.
% Connect MATLAB to a COMSOL server first; see docs/workflow.md.
root=setup_project;
if nargin<2,cfg=mechanical_config;end
if ischar(input)||isstring(input)
    data=load(input,'result');input=data.result;
end
pattern=prepare_mechanical_pattern(input);
if isempty(cfg.ligament_gap),cfg.ligament_gap=0.03+0.005*strcmp(pattern.mode,'rigid');end
if isempty(cfg.displacement),cfg.displacement=pattern.target.clampx-pattern.cols/2;end
assert(abs(cfg.displacement-(pattern.target.clampx-pattern.cols/2))<1e-8, ...
    'Prescribed grip displacement must match the target grip span.');
assert(exist('mphselectbox','file')~=0,'Add COMSOL LiveLink to MATLAB and connect with mphstart first.');
validateattributes(cfg.ga.radius,{'numeric'},{'scalar','finite','positive'});
validateattributes(cfg.ga.tolerance,{'numeric'},{'scalar','finite','nonnegative'});
old=rng;clean=onCleanup(@()rng(old));rng(cfg.seed,'twister');
out=tempname(fullfile(root,'results'));mkdir(out);
x0=pattern.x0;lower=x0-cfg.ga.radius;upper=x0+cfg.ga.radius;
if strcmp(pattern.mode,'rigid')
    lower(9:12)=max(0.02,lower(9:12));upper(9:12)=min(0.98,upper(9:12));
end
save(fullfile(out,'inputs.mat'),'pattern','cfg');
evaluate=@(x)evaluate_mechanical(x,pattern,cfg);
result=genetic_search(evaluate,x0,lower,upper,cfg.ga,@save_checkpoint);
result.pattern=pattern;result.config=cfg;result.output_dir=out;
save(fullfile(out,'mechanical_result.mat'),'result');
fprintf('Stopped: %s. Results: %s\n',result.stop_reason,out);
    function save_checkpoint(progress)
        save(fullfile(out,'checkpoint.mat'),'progress');
    end
end
