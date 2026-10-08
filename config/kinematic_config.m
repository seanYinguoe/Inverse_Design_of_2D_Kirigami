function cfg = kinematic_config
%KINEMATIC_CONFIG Shared settings for main_rigid and main_nonrigid.
cfg.rows = 2;
cfg.cols = 4;
cfg.target = make_target('vase','Size',3.2,'ClampX',2.5);
cfg.seed = 1;
cfg.make_plot = true;
% Fixed sheet size and symmetry match the original COMSOL parameterisation.
cfg.unit_length = 1;
cfg.solver = struct('FreeScale',false,'Symmetry',true,'EarlyStop',true, ...
    'Display','off','FeasibilityTolerance',1e-6,'Restarts',3);
end
