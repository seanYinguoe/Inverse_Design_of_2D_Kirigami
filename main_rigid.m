function result = main_rigid(cfg)
%MAIN_RIGID Kinematic inverse design; edit config/kinematic_config.m.
setup_project;
if nargin < 1, cfg = kinematic_config; end
result = run_kinematic_design('rigid',cfg);
end
