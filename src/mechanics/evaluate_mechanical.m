function [fitness,detail] = evaluate_mechanical(x,pattern,cfg)
%EVALUATE_MECHANICAL Rebuild, load and evaluate one candidate; release its model.
import com.comsol.model.util.*
nodes=decode_pattern(x,pattern);
if strcmp(pattern.mode,'rigid')
    [whole,half]=cuts_rigid(nodes,cfg.ligament_gap);
else
    [whole,half]=cuts_nonrigid(nodes,cfg.ligament_gap);
end
assert(isreal(whole)&&isreal(half)&&all(isfinite([whole(:);half(:)])), ...
    'Invalid cut geometry.');
model=build_comsol_model(whole,half,pattern,cfg);
cleanup=onCleanup(@()ModelUtil.remove(char(model.tag)));
model.sol('sol1').runAll;
s=cfg.metres_per_unit;e=1e-7*s;m=pattern.rows;n=pattern.cols;
box=[-n*s/2+e n*s/2-e;-m*s/2-e -m*s/2+e];
ids=mphselectbox(model,'geom1',box,'point');
assert(~isempty(ids),'No lower boundary points were found.');
% X,Y are reference coordinates, u,v are displacement: avoid frame ambiguity.
xx=mphevalpoint(model,'X+u','selection',ids,'solnum','end');
yy=mphevalpoint(model,'Y+v','selection',ids,'solnum','end');
boundary=[xx(:) yy(:)]/s;
[fitness,curve]=profile_error(boundary,pattern.target,cfg.fit_degree);
detail=struct('boundary',boundary,'sampled_x_range',[min(boundary(:,1)) max(boundary(:,1))], ...
    'fitted_boundary',curve,'nodes',{nodes});
end
