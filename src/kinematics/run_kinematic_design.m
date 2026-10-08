function result = run_kinematic_design(mode,cfg)
%RUN_KINEMATIC_DESIGN Shared runner; the rigid/nonrigid constraints stay separate.
root = setup_project;
assert(exist('fmincon','file') == 2,'Optimization Toolbox is required.');
validateattributes(cfg.rows,{'numeric'},{'scalar','integer','positive'});
validateattributes(cfg.cols,{'numeric'},{'scalar','integer','positive'});
if cfg.solver.Symmetry
    assert(all(mod([cfg.rows cfg.cols],2)==0),'Symmetric designs need even grid dimensions.');
end
cfg.mode = validatestring(mode,{'rigid','nonrigid'});
oldRng = rng; cleanup = onCleanup(@() rng(oldRng)); rng(cfg.seed,'twister');
[guess, initial] = fit_initial_guess(cfg.rows,cfg.cols,cfg.target,cfg.target.r,mode,cfg.unit_length);
opts = cfg.solver;
deployability = 1 + strcmp(mode,'nonrigid');
tic;
[deployed, diagnostics] = tessellation_optimization( ...
    guess,cfg.target,cfg.target.r,deployability,opts);
elapsed_seconds = toc;
compact = tessellation_compaction(deployed);
out = tempname(fullfile(root,'results')); mkdir(out);
result = struct('config',cfg,'solver_options',opts,'initial_guess',initial, ...
    'guess',{guess},'deployed',{deployed},'compact',{compact}, ...
    'diagnostics',diagnostics,'elapsed_seconds',elapsed_seconds, ...
    'matlab_version',version,'output_dir',out);
save(fullfile(out,'kinematic_result.mat'),'result');
if cfg.make_plot
    fig = figure('Color','w','Position',[100 100 1150 370]);
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    nexttile; plot_tessellation(guess); plot_boundary(cfg.target,cfg.target.r);
    axis off; title('Initial guess');
    nexttile; plot_tessellation(deployed); plot_boundary(cfg.target,cfg.target.r);
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
