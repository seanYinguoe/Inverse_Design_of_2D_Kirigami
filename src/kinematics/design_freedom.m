function [freedom, report] = design_freedom(m, n, s, r, p, opts)
%DESIGN_FREEDOM  How much design freedom is left after the constraints?
%
%   freedom = design_freedom(m, n, s, r, p)
%   [freedom, report] = design_freedom(m, n, s, r, p, opts)
%
%   RUN THIS BEFORE A LONG SOLVE. It answers, in a few seconds, the question
%   that otherwise costs an hour of fmincon: is this problem even solvable?
%
%   The design variables are the 2*m*n*16 node coordinates. Every equality
%   condition removes one dimension from the feasible set, PROVIDED it is
%   independent of the ones already imposed. This function builds the
%   Jacobian of all the equalities at a generic configuration and reports,
%   block by block, how many dimensions each one actually removes.
%
%   Read the result like this:
%     freedom > 0   the boundary conditions have room to be satisfied and
%                   the objective has something left to optimise over.
%     freedom == 0  the constraints pin the geometry completely. A solution
%                   exists only by coincidence; fmincon will grind and return
%                   a least-squares compromise that satisfies nothing exactly.
%                   The tell-tale sign is the BOUNDARY block adding 0 rank.
%     freedom < 0   cannot happen - rank is capped at the variable count -
%                   but a boundary block adding no rank means the same thing.
%
%   INPUTS
%     m, n : number of units along y and x
%     s    : target shape selector, see shape.m
%     r    : characteristic size of the target shape
%     p    : 1 = rigid deployable, 2 = non-rigid deployable
%     opts : same options struct as rigid.m / nonrigid.m (Symmetry,
%            FreeScale, MinEdgeLength). Use it to see what turning the
%            symmetry conditions on or off costs you.
%
%   OUTPUTS
%     freedom : dimensions left after all equalities (variables - rank)
%     report  : struct array, one entry per block, with name, rows, rank_added
%               and freedom_after
%
%   EXAMPLE
%     design_freedom(4,4,1,3.05,1)                       % rigid, circle
%     design_freedom(4,4,1,3.05,1,struct('Symmetry',true))
%
%   The Jacobian is built by finite differences over every variable, so cost
%   grows as (m*n)^2. It is fine up to about 6-by-6; beyond that, trust the
%   pattern rather than waiting.

if nargin < 6 || isempty(opts); opts = struct(); end
if ~isfield(opts,'Symmetry');      opts.Symmetry      = false; end
if ~isfield(opts,'FreeScale');     opts.FreeScale     = false; end
if ~isfield(opts,'MinEdgeLength'); opts.MinEdgeLength = 0;     end

if p == 1; fh = @rigid; mode = 'rigid'; else; fh = @nonrigid; mode = 'nonrigid'; end

N  = m*n*16;
nv = 2*N;

% a generic (slightly perturbed) configuration, so we measure the typical
% rank rather than one inflated or deflated by a special symmetry
rng(0);
T = tessellation_deployment(m, n, 1, pi/5);
x = units_to_nodes(T) + 1e-3*randn(N,2);

f = @(v) getceq(fh, v, s, m, n, r, opts);
ce0 = f(x(:));
J = zeros(numel(ce0), nv);
h = 1e-7;
for q = 1:nv
    xp = x(:); xp(q) = xp(q) + h;
    J(:,q) = (f(xp) - ce0)/h;
end

Aeq  = full(hinge_matrix(m,n));
acc  = Aeq;
rprev = rank(acc);

report = struct('name',{},'rows',{},'rank_added',{},'freedom_after',{});
report(end+1) = mk('hinges (linear, Aeq)', size(Aeq,1), rprev, nv-rprev);

B = block_map(m, n, mode, s, opts);
for q = 1:numel(B)
    idx = B{q}.idx;
    idx = idx(idx <= size(J,1));
    if isempty(idx); continue; end
    acc = [acc; J(idx,:)]; %#ok<AGROW>
    rnow = rank(acc);
    report(end+1) = mk(B{q}.name, numel(idx), rnow-rprev, nv-rnow); %#ok<AGROW>
    rprev = rnow;
end
freedom = nv - rprev;

fprintf('\n=== design_freedom : %s, shape %s, %d-by-%d, Symmetry=%d, FreeScale=%d ===\n', ...
    mode, describe_target(s), m, n, opts.Symmetry, opts.FreeScale);
fprintf('%d variables, %d equality rows\n\n', nv, size(acc,1));
fprintf('%-36s %6s %11s %14s\n','block','rows','rank added','freedom after');
for q = 1:numel(report)
    note = '';
    if report(q).rank_added < report(q).rows
        note = sprintf('  (%d redundant)', report(q).rows - report(q).rank_added);
    end
    if report(q).freedom_after == 0 && report(q).rank_added == 0
        note = [note '  <-- CANNOT BE SATISFIED'];
    end
    fprintf('%-36s %6d %11d %14d%s\n', report(q).name, report(q).rows, ...
        report(q).rank_added, report(q).freedom_after, note);
end
fprintf('\nremaining design freedom: %d\n', freedom);
if freedom <= 0
    fprintf(['VERDICT: over-constrained. The optimiser cannot satisfy the target\n' ...
             '         boundary. Turn Symmetry off, use FreeScale, or use a\n' ...
             '         larger grid to buy back degrees of freedom.\n']);
else
    fprintf('VERDICT: solvable - %d dimensions left for the objective.\n', freedom);
end
end

% -------------------------------------------------------------------------
function str = describe_target(s)
% a printable name for either a built-in code or a make_target struct
if isstruct(s)
    switch s.type
        case 'implicit'; str = 'custom (implicit)';
        case 'curve';    str = sprintf('custom (%d-point curve)', size(s.pts,1));
        case 'builtin';  str = sprintf('%d', s.code);
        otherwise;       str = 'custom';
    end
    if isfield(s,'clampx') && ~isempty(s.clampx)
        str = [str sprintf(', clamped at x=+-%g', s.clampx)];
    end
else
    str = sprintf('%d', s);
end
end

% -------------------------------------------------------------------------
function e = mk(name, rows, added, freedom)
e = struct('name',name,'rows',rows,'rank_added',added,'freedom_after',freedom);
end

function v = getceq(fh, x, s, m, n, r, opts)
[~, ceq] = fh(reshape(x, [], 2), s, m, n, r, opts);
v = ceq(:);
end

function Aeq = hinge_matrix(m,n)
N = m*n*16; nv = 2*N;
index1 = [4 10 6 16 8 3 15 2];
index2 = [7 13 11 1 9 14 12 5];
rows = []; cols = []; vals = []; e = 1;
    function push(a,b)
        rows(end+1:end+4) = [e e e+1 e+1];
        cols(end+1:end+4) = [a b a+N b+N];
        vals(end+1:end+4) = [1 -1 1 -1];
        e = e + 2;
    end
for i=1:m, for j=1:n, for k=1:4
    push(index1(k)+(i-1)*n*16+16*(j-1), index2(k)+(i-1)*n*16+16*(j-1)); end, end, end
for i=1:m, for j=1:n-1, for k=5:6
    push(index1(k)+(i-1)*n*16+16*(j-1), index2(k)+(i-1)*n*16+16*j); end, end, end
for i=1:m-1, for j=1:n, for k=7:8
    push(index1(k)+(i-1)*n*16+16*(j-1), index2(k)+i*n*16+16*(j-1)); end, end, end
Aeq = sparse(rows, cols, vals, e-1, nv);
end

function B = block_map(m, n, mode, ~, opts)
% index ranges of each ceq block, mirroring the loop structure of
% rigid.m / nonrigid.m. Keep in step with those files.
L = {};
A = @(name,cnt) struct('name',name,'cnt',cnt);
if strcmp(mode,'rigid')
    L{end+1} = A('angle in-unit (2pi + rigid pi)', m*n*3);
    L{end+1} = A('angle rigid pi across rows',     (m-1)*n);
else
    L{end+1} = A('angle in-unit 2pi',              m*n);
end
L{end+1} = A('angle adj horiz 2pi', m*(n-1));
L{end+1} = A('angle adj vert 2pi',  n*(m-1));
L{end+1} = A('angle 4-unit corner', (m-1)*(n-1));
if strcmp(mode,'rigid')
    L{end+1} = A('edge in-unit (parallelogram)', m*n*2);
    L{end+1} = A('edge adj horiz',        m*(n-1)*2);
    L{end+1} = A('edge adj vert',         (m-1)*n*2);
    L{end+1} = A('edge diagonal',         (m-1)*(n-1)*2);
    L{end+1} = A('edge border rows/cols',  n+n+(m-1)+(m-1));
else
    L{end+1} = A('edge in-unit',   m*n*4);
    L{end+1} = A('edge adj horiz', m*(n-1)*2);
    L{end+1} = A('edge adj vert',  (m-1)*n*2);
end
if opts.Symmetry
    L{end+1} = A('symmetry left-right', (m/2)*(n/2)*16*2);
    L{end+1} = A('symmetry top-bottom', (m/2)*(n/2)*8*2);
end
if opts.FreeScale
    L{end+1} = A('square size (aspect only)', 3);
else
    L{end+1} = A('square size + aspect', 4);
end
L{end+1} = A('square edge straightness', n+n+m+m+(m-1)+(m-1)+(n-1)+(n-1));
L{end+1} = A('square corner right angles', 4);
L{end+1} = A('BOUNDARY target curve', 4*m + 4*n);

B = {}; k = 1;
for q = 1:numel(L)
    B{q}.name = L{q}.name;
    B{q}.idx  = k:(k+L{q}.cnt-1);
    k = k + L{q}.cnt;
end
end
