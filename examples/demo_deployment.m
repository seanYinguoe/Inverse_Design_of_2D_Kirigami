function result = demo_deployment(makePlot)
%DEMO_DEPLOYMENT Rotate a 4-by-4 sheet; no optimisation toolbox is required.
% The animation is kinematic: it does not model elastic hinge forces.
if nargin < 1, makePlot = true; end
root = setup_project;
angles = [0, pi/12, pi/4];
states = cell(size(angles));
for k = 1:numel(angles)
    states{k} = tessellation_deployment(4, 4, 1, angles(k));
end
out = tempname(fullfile(root, 'results')); mkdir(out);
result = struct('angles', angles, 'states', {states}, 'output_dir', out);
save(fullfile(out, 'deployment.mat'), 'result');
if makePlot
    fig = figure('Color','w','Position',[100 100 1150 360]);
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    for k = 1:numel(angles)
        nexttile; plot_tessellation(states{k}); axis off;
        title(sprintf('Panel rotation: %.0f degrees', angles(k)*180/pi));
    end
    exportgraphics(fig, fullfile(out,'deployment.png'), 'Resolution',180);
end
fprintf('Deployment example saved to %s\n', out);
end
