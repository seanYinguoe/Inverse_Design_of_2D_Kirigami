function result = main_nonrigid(cfg)
%MAIN_NONRIGID Kinematic inverse design; edit config/kinematic_config.m.
setup_project;
if nargin < 1, cfg = kinematic_config; end
result = run_kinematic_design('nonrigid',cfg);
end
