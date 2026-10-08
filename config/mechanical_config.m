function cfg = mechanical_config
%MECHANICAL_CONFIG Explicit example settings inherited from the source scripts.
% Review units and material values for your specimen; these are not a calibration.
cfg.material = struct('density',2700,'young_modulus',4.33e9,'poisson_ratio',0.2);
cfg.metres_per_unit = 1;     % original model used metres; set your physical scale
cfg.thickness = 0.003;      % metres
cfg.cut_width = 0.02;       % design-coordinate units
cfg.ligament_gap = [];      % source defaults: rigid .035, nonrigid .03
cfg.displacement = [];      % per grip, derived from target.clampx - cols/2
cfg.load_steps = 50;
cfg.mesh_size = 6;          % COMSOL automatic size: 1 finest, 9 coarsest
cfg.fit_degree = 5;
cfg.seed = 1;
cfg.ga = struct('population',20,'generations',15,'crossover',0.6, ...
    'mutation',0.01,'radius',0.08,'tolerance',0.005,'stall_generations',5);
end
