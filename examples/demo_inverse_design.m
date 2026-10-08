function result = demo_inverse_design(mode, makePlot)
%DEMO_INVERSE_DESIGN Stage-1 circle example, with explicit solver diagnostics.
% result = demo_inverse_design('rigid') or demo_inverse_design('nonrigid').
% This example uses the current maintained solver, not a frozen paper run.
if nargin < 1, mode = 'rigid'; end
if nargin < 2, makePlot = true; end
mode = validatestring(mode, {'rigid','nonrigid'});
root = setup_project;
assert(exist('fmincon','file') == 2, 'Optimization Toolbox is required.');
cfg = struct('rows',4,'cols',4,'target_shape',1, ...
    'shape_size',3.2*sqrt(2/2.2),'mode',mode,'seed',1);
oldRng = rng; cleanup = onCleanup(@() rng(oldRng)); rng(cfg.seed,'twister');
[guess, initial] = fit_initial_guess(cfg.rows,cfg.cols,cfg.target_shape,cfg.shape_size,mode);
opts = struct('FreeScale',true,'EarlyStop',true,'Display','off', ...
    'FeasibilityTolerance',1e-4,'Restarts',3);
deployability = 1 + strcmp(mode,'nonrigid');
tic;
[deployed, diagnostics] = tessellation_optimization( ...
    guess,cfg.target_shape,cfg.shape_size,deployability,opts);
elapsed_seconds = toc;
compact = tessellation_compaction(deployed);
out = tempname(fullfile(root,'results')); mkdir(out);
result = struct('config',cfg,'solver_options',opts,'initial_guess',initial, ...
    'guess',{guess},'deployed',{deployed},'compact',{compact}, ...
    'diagnostics',diagnostics,'elapsed_seconds',elapsed_seconds, ...
    'matlab_version',version,'output_dir',out);
save(fullfile(out,'inverse_design.mat'),'result');
if makePlot
    fig = figure('Color','w','Position',[100 100 1150 370]);
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    nexttile; plot_tessellation(guess); plot_boundary(cfg.target_shape,cfg.shape_size);
    axis off; title('Initial guess');
    nexttile; plot_tessellation(deployed); plot_boundary(cfg.target_shape,cfg.shape_size);
    axis off; title(sprintf('Designed boundary (residual %.1e)',diagnostics.maxceq));
    nexttile; plot_tessellation(compact); axis off; title('Compacted geometry');
    exportgraphics(fig,fullfile(out,'inverse_design.png'),'Resolution',180);
end
fprintf('%s: converged=%d, max equality=%.3g, inequality=%.3g, hinge=%.3g; %.1f s\n', ...
    mode,diagnostics.converged,diagnostics.maxceq,diagnostics.maxc,diagnostics.maxlin,elapsed_seconds);
fprintf('Saved configuration and diagnostics to %s\n',out);
if ~diagnostics.converged
    warning('kirigami:InfeasibleDesign','Design is not feasible at the stated tolerance. Inspect diagnostics before using it.');
end
end
