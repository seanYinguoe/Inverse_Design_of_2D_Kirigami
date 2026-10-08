function [tessellation_optimized, info] = tessellation_optimization(tessellation_transformed,s,r,p,opts)
%TESSELLATION_OPTIMIZATION  Stage-1 kinematic optimisation of the cut pattern.
%
%   tessellation_optimized = tessellation_optimization(guess,s,r,p)
%   [tessellation_optimized, info] = tessellation_optimization(guess,s,r,p,opts)
%
%   Moves every node of the tessellation with fmincon until the deployed
%   boundary sits on the target shape, subject to the geometric conditions of
%   Section 2.1 of the paper. The design variables are the coordinates of all
%   m*n*16 nodes; the hinges are imposed as LINEAR equalities (Aeq) and the
%   deployability / non-overlap / boundary conditions as NONLINEAR ones
%   (rigid.m or nonrigid.m).
%
%   INPUTS
%     tessellation_transformed : m-by-n cell array, the deployed pattern used
%                                as the initial guess. Use fit_initial_guess
%                                rather than a hand-picked opening angle.
%     s : target shape selector, see shape.m (1 circle ... 5 heart)
%     r : characteristic size of the target shape
%     p : deployability mode, 1 = rigid deployable, 2 = non-rigid deployable
%     opts : optional struct, see below
%
%   OPTIONS (all optional)
%     MaxFunctionEvaluations : default 400*nvars. THIS IS THE ONE THAT
%           MATTERS. fmincon builds gradients by finite differences, so one
%           iteration costs nvars+1 evaluations. The old hard-coded 8000 gave
%           a 4x4 grid (512 variables) about FIFTEEN iterations before the
%           solver gave up - nowhere near enough to satisfy 465 nonlinear
%           equalities, which is why the compacted state came out warped.
%     MaxIterations          : default 3000
%     Restarts               : default 3. If fmincon stops on an iteration or
%           evaluation limit, it is restarted from where it stopped. A warm
%           restart resets the barrier parameter and often clears a stall.
%     MinEdgeLength          : default 0.15*unit size. Floor on panel edges,
%           see rigid.m. Set 0 to disable.
%     Symmetry               : default FALSE. Mirror-symmetry conditions.
%           They are off by default because they over-constrain the problem:
%           for a 4x4 rigid circle they consume so much of the design freedom
%           that the target-boundary equations add no rank at all and simply
%           cannot be met. Check the returned feasibility residuals for your case.
%     FreeScale              : default false. Let the compacted sheet size
%           float, pinning only its aspect ratio, so the optimiser can scale
%           the design to fit the boundary. Pair with fit_initial_guess.
%     FeasibilityTolerance   : default 1e-7. How small the constraint
%           violations must be to call the result converged, and what is
%           handed to fmincon as its ConstraintTolerance. Do not set this
%           much tighter: on a unit-scale geometry 1e-7 is already far below
%           any manufacturing tolerance, and asking for 1e-8 makes fmincon
%           report exitflag -2 ("no feasible point") on solutions whose worst
%           violation is ~1e-7 and which are perfectly good designs.
%     UseParallel            : default FALSE. fmincon can spread its
%           finite-difference evaluations over a parallel pool, and since one
%           gradient costs nvars+1 full evaluations that looks like an easy
%           win. MEASURED, IT IS NOT: on a 4x4 rigid circle with 15 cores it
%           came out at 0.6x, i.e. 1.6x SLOWER than serial. One constraint
%           evaluation only takes ~16 ms, which is too little work to cover
%           parfor's per-iteration dispatch overhead. Try it on much larger
%           grids, where each evaluation is dearer, but measure before
%           trusting it.
%     SparseGradients        : default FALSE. Build the objective and
%           constraint Jacobians by SPARSE finite differences (sparse_fd.m)
%           instead of letting fmincon perturb every variable.
%
%           The sparsity argument is sound and the saving is real: each
%           constraint touches only ~6 of the variables, so 73 colours suffice
%           at 4x4 where fmincon uses 513 evaluations, and the ratio improves
%           with grid size. Measured on a 4x4 rigid circle the evaluation
%           count fell from 186755 to between 1765 and 9335.
%
%           IT DOES NOT PAY OFF, because evaluations are not the whole cost.
%           Supplying gradients changes fmincon's search path, and the
%           iteration count swings enough to swamp the saving:
%
%             dense (fmincon)  186755 evals   360 iters   191 s  converged
%             sparse forward     1765 evals   318 iters    94 s  NOT converged
%             sparse central     9335 evals  3104 iters  1024 s  converged
%
%           Forward differences are twice as fast but stall at max|ceq| ~2e-7,
%           above the feasibility tolerance. Central differences converge but
%           take five times as long. Left off by default; switch it on only if
%           you are prepared to measure your own case.
%
%           One genuine side effect worth knowing: the sparse-gradient runs
%           consistently found BETTER optima (objective 0.649 and 0.682 versus
%           0.816 for the dense path). If solution quality matters more than
%           run time, it is worth trying.
%     ConstraintTolerance    : default 1e-6. How feasible a point must be
%           before fmincon will call it acceptable. NOTE: on this problem it
%           does NOT control run time - the solver exits on StepTolerance
%           (exitflag 2, "step size below tolerance, constraints satisfied"),
%           so it never gets as far as consulting this. Sweeping it from 1e-5
%           to 1e-8 changed nothing: same 360 iterations, same answer.
%           1e-6 is chosen because it is physically ample - the model is ~4
%           units across, so 1e-6 is a nanometre on a 100 mm sheet.
%     OptimalityTolerance    : default 1e-6, same story.
%     FeasibilityTolerance   : default 1e-6. The threshold info.converged is
%           judged against. Keep it at or above ConstraintTolerance.
%     StepTolerance          : default [] (fmincon's own, 1e-10). Without
%           EarlyStop this is the criterion that actually ends the solve.
%     EarlyStop              : default TRUE. Stop as soon as the design is
%           FEASIBLE and the objective has stopped meaningfully improving,
%           rather than grinding to fmincon's step tolerance. On a 4x4 circle
%           this is the difference between ~20 s and ~190 s.
%
%           A fixed MaxIterations cap does not work well as a substitute:
%           interior-point feasibility is not monotone, so a cap can land on
%           an iterate that happens to be infeasible (capping at 100 gave
%           |ceq| = 2.6e-05, worse than capping at 60). EarlyStop only ever
%           stops at a point it has checked is feasible.
%     EarlyStopWindow        : default 15. Iterations of objective history
%           used to judge "stopped improving".
%     EarlyStopRelImprove    : default 0.01. Stop when the objective has
%           improved by less than this fraction over the window.
%     MaxIterations          : default 3000, and the blunt way to cap the run.
%
%           WHERE THE TIME GOES. On a 4x4 circle the geometry is essentially
%           solved early and the rest is objective polishing:
%
%               |ceq| < 1e-04   iteration  25    13 s
%               |ceq| < 1e-06   iteration  33    16 s
%               |ceq| < 1e-07   iteration 265   140 s
%               stops                     360   195 s
%
%           So a feasible design exists after ~16 s. Everything after that is
%           making the panels more uniform, which is what the objective asks
%           for. Capping MaxIterations trades run time against that uniformity
%           - it does not make the design invalid. Check info.maxceq to see
%           where you landed.
%     Display                : fmincon display, default 'iter'.
%
%   OUTPUTS
%     tessellation_optimized : m-by-n cell array of 16-by-2 node lists
%     info : struct with exitflag, fval, iterations, funcCount, attempts,
%            maxceq / maxc (worst constraint violation at the solution),
%            converged (logical), and the achieved compact size.
%
%   ALWAYS CHECK info.converged. An unconverged result still looks like a
%   tessellation but does not satisfy the geometric conditions, so
%   tessellation_compaction will not fold it into a proper rectangle.
%
%   Feed the result to tessellation_compaction to obtain the as-cut pattern.

[m , n] = size(tessellation_transformed);
N = m*n*16;
nvars = 2*N;

%% Options
if nargin < 5 || isempty(opts); opts = struct(); end
unit_scale = estimate_unit_size(tessellation_transformed);
if ~isfield(opts,'MaxFunctionEvaluations'); opts.MaxFunctionEvaluations = 400*nvars; end
if ~isfield(opts,'MaxIterations');          opts.MaxIterations = 3000;               end
if ~isfield(opts,'Restarts');               opts.Restarts = 3;                       end
opts.Restarts = max(0, round(opts.Restarts));   % the solve loop must run at least once
if ~isfield(opts,'MinEdgeLength');          opts.MinEdgeLength = 0.15*unit_scale;    end
if ~isfield(opts,'FreeScale');              opts.FreeScale = false;                  end
if ~isfield(opts,'Symmetry');               opts.Symmetry = false;                   end
if ~isfield(opts,'FeasibilityTolerance');   opts.FeasibilityTolerance = 1e-6;        end
if ~isfield(opts,'UseParallel');            opts.UseParallel = false;                end
if ~isfield(opts,'SparseGradients');        opts.SparseGradients = false;            end
if ~isfield(opts,'ConstraintTolerance');    opts.ConstraintTolerance = 1e-6;         end
if ~isfield(opts,'OptimalityTolerance');    opts.OptimalityTolerance = 1e-6;         end
if ~isfield(opts,'StepTolerance');          opts.StepTolerance = [];                 end
if ~isfield(opts,'EarlyStop');              opts.EarlyStop = true;                   end
if ~isfield(opts,'EarlyStopWindow');        opts.EarlyStopWindow = 15;               end
if ~isfield(opts,'EarlyStopRelImprove');    opts.EarlyStopRelImprove = 0.01;         end
if ~isfield(opts,'Display');                opts.Display = 'iter';                   end
if ~isfield(opts,'OutputFcn');              opts.OutputFcn = [];                     end

conopts = struct('MinEdgeLength',opts.MinEdgeLength, ...
                 'FreeScale',opts.FreeScale, ...
                 'Symmetry',opts.Symmetry);

%% Define the objective function
unit_length = 1;
% number of interior corners the two objective terms are summed over,
% used to normalise the objective so it does not scale with the grid size
M = (m-1)*(n-1)*2 + (m-1)+(n-1);
% keep the pattern as close to periodic as the boundary condition allows:
% angle_diff has units of rad^2 and is converted to a length by unit_length/pi
% Written as a residual vector r with f = sum(r.^2), so that the gradient is
% 2*J'*r and J can be built sparsely (see objective_residuals.m).
w_angle = sqrt(unit_length/(M*pi));
w_edge  = sqrt(1/M);
resfun  = @(xv) objective_residuals(xv, m, n, w_angle, w_edge);
fun     = @(nodes) 1/M * (unit_length/pi * (angle_diff(nodes,m,n)) + edge_diff(nodes,m,n));

%% Define the initial configuration
% fmincon optimises one flat (m*n*16)-by-2 array of node coordinates
tessellation_initial = units_to_nodes(tessellation_transformed);

%% Define the nonlinear conditions (deployability, non-overlap, boundary)
if p == 1
    nonlcon = @(nodes) rigid(nodes,s,m,n,r,conopts);
elseif p == 2
    nonlcon = @(nodes) nonrigid(nodes,s,m,n,r,conopts);
else
    error('tessellation_optimization:badMode','p must be 1 (rigid) or 2 (non-rigid).');
end
A = []; % linear inequality constraints
b = []; % linear inequality constraints
lb = []; % Lower bounds
ub = []; % Upper bounds

%% Hinge conditions as linear equalities: Aeq*x = 0
% fmincon flattens the (m*n*16)-by-2 variable column-major, so the unknown
% vector is [all x coordinates; all y coordinates] and the y coordinate of
% node k lives at column k + m*n*16. Each hinge contributes two rows
% (x and y): "these two duplicated nodes must stay at the same place".
%
% Aeq is built SPARSE and then trimmed to the rows actually used. The old
% version allocated a dense nvars-by-nvars matrix and filled only a fraction
% of it - for a 4x4 grid that is 224 real rows and 1824 all-zero ones. Those
% zero rows are trivially satisfied but leave Aeq badly rank deficient, and
% the interior-point KKT solve inherits that (RCOND ~ 1e-19 in testing).
% It is also O(grid^4) memory: 5x5 is a 800x800 dense matrix, 10x10 is
% 3200x3200.
% hinge node pairs (see the node map in create_unit.m)
%   k = 1..4 : hinges inside a unit   4-7 right-mid, 10-13 left-mid,
%                                     6-11 and 16-1 at the centre
%   k = 5..6 : hinges to the unit on the right  8-9 and 3-14
%   k = 7..8 : hinges to the unit above        15-12 and 2-5
index1 = [4 10 6 16 8 3 15 2];
index2 = [7 13 11 1 9 14 12 5];

nHinge = (m*n*4 + m*(n-1)*2 + (m-1)*n*2);   % hinges, each giving an x and a y row
rows = zeros(1, 4*nHinge);
cols = zeros(1, 4*nHinge);
vals = zeros(1, 4*nHinge);
e = 0;      % running index into the triplet lists
eq_row = 1; % running row of Aeq

    function push(a_col, b_col)
        % one hinge, x row then y row
        rows(e+1) = eq_row;   cols(e+1) = a_col;     vals(e+1) =  1;
        rows(e+2) = eq_row;   cols(e+2) = b_col;     vals(e+2) = -1;
        rows(e+3) = eq_row+1; cols(e+3) = a_col + N; vals(e+3) =  1;
        rows(e+4) = eq_row+1; cols(e+4) = b_col + N; vals(e+4) = -1;
        e = e + 4;
        eq_row = eq_row + 2;
    end

% hinges internal to each unit
for i = 1:m
    for j = 1:n
        for k = 1:4
            push(index1(k)+(i-1)*n*16+16*(j-1), index2(k)+(i-1)*n*16+16*(j-1));
        end
    end
end
% hinges between horizontally adjacent units (j and j+1)
for i = 1:m
    for j = 1:n-1
        for k = 5:6
            push(index1(k)+(i-1)*n*16+16*(j-1), index2(k)+(i-1)*n*16+16*j);
        end
    end
end
% hinges between vertically adjacent units (i and i+1)
for i = 1:m-1
    for j = 1:n
        for k = 7:8
            push(index1(k)+(i-1)*n*16+16*(j-1), index2(k)+i*n*16+16*(j-1));
        end
    end
end
Aeq = sparse(rows, cols, vals, eq_row-1, nvars);
beq = zeros(eq_row-1, 1);

%% Sparse gradient machinery
if opts.SparseGradients
    x0v = tessellation_initial(:);
    confun_vec = @(xv) con_vector(nonlcon, xv);
    if ~strcmp(opts.Display,'off')
        fprintf('analysing sparsity...\n');
    end
    Pc = sparse_fd('prepare', confun_vec, x0v);
    Pr = sparse_fd('prepare', resfun,     x0v);
    [nc_, ~] = nonlcon(tessellation_initial);
    nIneq = numel(nc_);
    if ~strcmp(opts.Display,'off')
        fprintf(['  constraints: %d rows, %.2f%% dense, %d colours ' ...
                 '(vs %d dense evaluations, %.1fx)\n'], ...
            Pc.nrows, Pc.density*100, Pc.ncolour, nvars+1, (nvars+1)/(Pc.ncolour+1));
        fprintf(['  objective  : %d residuals, %.2f%% dense, %d colours ' ...
                 '(vs %d dense evaluations, %.1fx)\n'], ...
            Pr.nrows, Pr.density*100, Pr.ncolour, nvars+1, (nvars+1)/(Pr.ncolour+1));
    end
    objective = @(xm) obj_with_grad(xm, resfun, Pr, M);
    constraints = @(xm) con_with_grad(xm, nonlcon, confun_vec, Pc, nIneq);
else
    objective = fun;
    constraints = nonlcon;
end

%% Solve
% Early stop: leave as soon as the point is feasible AND the objective has
% flattened out. History is reset before each fmincon attempt.
histF = [];
    function stop = earlystop_fcn(~, ov, state)
        stop = false;
        if ~isempty(opts.OutputFcn)
            stop = opts.OutputFcn([], ov, state);
        end
        if strcmp(state,'init'); histF = []; return; end
        if ~strcmp(state,'iter') || ~opts.EarlyStop; return; end
        histF(end+1) = ov.fval;
        W = opts.EarlyStopWindow;
        if numel(histF) <= W; return; end
        if ov.constrviolation > opts.FeasibilityTolerance; return; end
        improved = (histF(end-W) - histF(end)) / max(1, abs(histF(end)));
        if improved < opts.EarlyStopRelImprove
            stop = true;
        end
    end
if opts.EarlyStop || ~isempty(opts.OutputFcn)
    outfcn = @earlystop_fcn;
else
    outfcn = [];
end

options = optimoptions('fmincon', ...
    'Display', opts.Display, ...
    'Algorithm', 'interior-point', ...
    'MaxFunctionEvaluations', opts.MaxFunctionEvaluations, ...
    'MaxIterations', opts.MaxIterations, ...
    'OptimalityTolerance', opts.OptimalityTolerance, ...
    'ConstraintTolerance', opts.ConstraintTolerance, ...
    'ScaleProblem', true, ...
    'SpecifyObjectiveGradient', opts.SparseGradients, ...
    'SpecifyConstraintGradient', opts.SparseGradients, ...
    'UseParallel', opts.UseParallel, ...
    'OutputFcn', outfcn);
if ~isempty(opts.StepTolerance)
    options = optimoptions(options, 'StepTolerance', opts.StepTolerance);
end

if opts.UseParallel && isempty(gcp('nocreate'))
    fprintf('starting a parallel pool for the finite-difference gradients...\n');
    parpool('threads');   % threads pool: no data copying, starts in ~1 s
end
if opts.UseParallel && ~has_parallel()
    warning('tessellation_optimization:noParallel', ...
        'UseParallel requested but the Parallel Computing Toolbox is unavailable.');
end

x = tessellation_initial;
totalFunc = 0; totalIter = 0;
for attempt = 1:(1 + opts.Restarts)
    [x, fval, exitflag, output] = fmincon(objective, x, A, b, Aeq, beq, lb, ub, constraints, options);
    totalFunc = totalFunc + output.funcCount;
    totalIter = totalIter + output.iterations;
    if exitflag > 0
        break            % converged
    end
    if opts.EarlyStop && exitflag == -1
        exitflag = 2;    % stopped by our OutputFcn at a feasible point
        break
    end
    if exitflag < 0
        break            % infeasible or failed - restarting will not help
    end
    % exitflag == 0: hit an iteration or evaluation limit, warm restart
end

% back from the flat node array to the m-by-n cell array of units
tessellation_optimized = nodes_to_units(x,m,n);

%% Report
[cf, ceqf] = nonlcon(x);
info.exitflag  = exitflag;
info.fval      = fval;
info.iterations = totalIter;
info.funcCount = totalFunc;
info.attempts  = attempt;
info.maxceq    = max(abs([ceqf(:); 0]));
info.maxc      = max([cf(:); 0]);
info.maxlin    = max(abs([Aeq*x(:) - beq; 0]));
% Judge on the violations actually achieved rather than on exitflag alone.
% fmincon returns -2 whenever it cannot reach ConstraintTolerance, even when
% the point it stopped at is a good design; conversely a positive flag with a
% large residual is not a usable result. The geometry is what matters.
tol = opts.FeasibilityTolerance;
info.tolerance = tol;
info.converged = info.maxceq <= tol && info.maxc <= tol && info.maxlin <= tol;

if ~strcmp(opts.Display,'off')
    fprintf('\n--- tessellation_optimization ---\n');
    fprintf('  exitflag %d after %d iteration(s) in %d attempt(s), %d evaluations\n', ...
        info.exitflag, info.iterations, info.attempts, info.funcCount);
    fprintf('  objective        %.6g\n', info.fval);
    fprintf('  max |ceq|        %.3e   (nonlinear equalities)\n', info.maxceq);
    fprintf('  max c            %.3e   (nonlinear inequalities, want <= 0)\n', info.maxc);
    fprintf('  max |Aeq x|      %.3e   (hinges)\n', info.maxlin);
    fprintf('  tolerance        %.3e\n', tol);
    if info.converged
        fprintf('  CONVERGED - the geometric conditions are satisfied.\n');
    else
        fprintf(['  NOT CONVERGED: inspect residuals, target and initial guess.\n' ...
                 '  The current design must not be sent to the mechanical stage.\n']);
    end
end
end

% -------------------------------------------------------------------------
function v = con_vector(nonlcon, xv)
% both constraint vectors stacked, as a column, for Jacobian probing
[c, ceq] = nonlcon(reshape(xv, [], 2));
v = [c(:); ceq(:)];
end

% -------------------------------------------------------------------------
function [f, g] = obj_with_grad(xm, resfun, Pr, ~)
% f = sum(r.^2), so df/dx = 2*J'*r with J the (sparse) Jacobian of r.
% The 1/M normalisation is already folded into the residual weights,
% so no further scaling is applied here.
r = resfun(xm(:));
f = sum(r.^2);
if nargout > 1
    J = sparse_fd('jacobian', Pr, resfun, xm(:));
    g = reshape(full(2*(J.'*r)), size(xm));
end
end

% -------------------------------------------------------------------------
function [c, ceq, gc, gceq] = con_with_grad(xm, nonlcon, confun_vec, Pc, nIneq)
% fmincon wants gradients as columns: gc is nvars-by-(number of constraints)
[c, ceq] = nonlcon(xm);
if nargout > 2
    J = sparse_fd('jacobian', Pc, confun_vec, xm(:));
    Jt = J.';                       % nvars-by-nrows
    gc   = full(Jt(:, 1:nIneq));
    gceq = full(Jt(:, nIneq+1:end));
end
end

% -------------------------------------------------------------------------
function tf = has_parallel()
% true when the Parallel Computing Toolbox is installed and licensed
tf = ~isempty(ver('parallel')) && license('test','Distrib_Computing_Toolbox');
end

% -------------------------------------------------------------------------
function L = estimate_unit_size(T)
% typical quadrant-square edge length of the guess, used to scale defaults
[m, n] = size(T);
P = T{max(1,round(m/2)), max(1,round(n/2))};
L = mean([norm(P(2,:)-P(1,:)), norm(P(3,:)-P(2,:)), ...
          norm(P(4,:)-P(3,:)), norm(P(1,:)-P(4,:))]);
if ~isfinite(L) || L <= 0; L = 0.5; end
end
