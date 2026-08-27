function out = sparse_fd(action, varargin)
%SPARSE_FD  Sparse finite-difference Jacobians by graph colouring.
%
%   Why this exists
%     fmincon has no analytic gradients for this problem, so it perturbs every
%     variable one at a time: nvars+1 evaluations for ONE Jacobian. For a 4x4
%     grid that is 513 evaluations of the whole constraint set per iteration.
%
%     But the constraint Jacobian is very sparse - each condition involves
%     about three nodes, so a row has ~6 nonzeros out of 512. Two variables
%     that never appear in the same constraint can be perturbed TOGETHER in a
%     single evaluation, and their derivatives read off independently, because
%     no row receives a contribution from both. The smallest number of such
%     groups is a colouring of the column-intersection graph.
%
%     Measured on this problem: 73 colours instead of 513 evaluations at 4x4,
%     109 instead of 1153 at 6x6. The denser the grid, the bigger the win.
%
%   USAGE
%     P = sparse_fd('prepare', f, x0)          analyse f at x0: probe the
%                                              sparsity pattern and colour it
%     J = sparse_fd('jacobian', P, f, x)       Jacobian of f at x, using P
%
%   f must take a column vector x and return a column vector of values whose
%   sparsity structure does not change with x.
%
%   ASSEMBLY  The Jacobian is assembled straight into triplet form: the
%   (row, column) lists are computed ONCE in 'prepare', and each build only
%   fills a values vector. An earlier version did a find() per column and
%   allocated a dense nrows-by-nvars matrix on every build, which cost more
%   than the evaluations it was saving.
%
%   STEP SIZE  Relative, h_j = h0 * max(1, |x_j|), matching what fmincon does
%   internally. A fixed absolute step gives noticeably worse gradients on
%   coordinates of varying magnitude, and a bad gradient costs more iterations
%   than it saves in evaluations.
%
%   FORWARD vs CENTRAL  Pass 'Type','central' to prepare (the default) for
%   two-sided differences: 2*ncolours evaluations instead of ncolours+1, with
%   error O(h^2) instead of O(h). That accuracy is not optional here. Measured
%   on a 4x4 rigid circle, one-sided differences ran 2.0x faster than dense
%   fmincon but stalled at max|ceq| ~ 2e-7, above the feasibility tolerance,
%   so the result never certified as converged. Two-sided costs half the
%   speed-up and reaches the same solution quality as the dense path.
%
%   THE PATTERN IS CACHED. It is probed once, at a perturbed configuration so
%   that no entry is accidentally zero at a symmetric point. If you change the
%   constraint set (different options, different shape) you must prepare again
%   - tessellation_optimization does this on every call.

switch lower(action)
    case 'prepare'
        out = do_prepare(varargin{:});
    case 'jacobian'
        out = do_jacobian(varargin{:});
    otherwise
        error('sparse_fd:badAction','Unknown action ''%s''.', action);
end
end

% -------------------------------------------------------------------------
function P = do_prepare(f, x0, varargin)
ip = inputParser;
ip.addParameter('Step', eps^(1/3), @(v) isnumeric(v) && isscalar(v));
ip.addParameter('Type', 'central', @(v) any(strcmpi(v,{'forward','central'})));
ip.parse(varargin{:});
h0   = ip.Results.Step;
ftype = lower(ip.Results.Type);
x0 = x0(:);
nv = numel(x0);

% Probe the pattern away from any symmetric configuration: at a symmetric
% point a derivative can vanish by cancellation and be mistaken for a
% structural zero, which would silently truncate the Jacobian.
rng(0);
xp = x0 + 1e-3*randn(nv,1);
v0 = f(xp);
nr = numel(v0);

hp = 1e-7;
rows = cell(1,nv);
cols = cell(1,nv);
for q = 1:nv
    xq = xp; xq(q) = xq(q) + hp;
    d = (f(xq) - v0)/hp;
    idx = find(abs(d) > 1e-8);
    rows{q} = idx(:).';
    cols{q} = repmat(q, 1, numel(idx));
end
S = sparse([rows{:}], [cols{:}], true, nr, nv);

colour = colour_columns(S);
ncol   = max(colour);

% triplet structure, fixed for the life of P
[Ir, Jc] = find(S);
Ir = Ir(:); Jc = Jc(:);

% which nonzeros belong to each colour group, and which columns each group
% perturbs - both precomputed so a build is pure vectorised assignment
groupCols = accumarray(colour(:), (1:nv).', [], @(v){v(:)});
nzColour  = colour(Jc).';
groupNz   = accumarray(nzColour(:), (1:numel(Ir)).', [], @(v){v(:)});

P.type      = ftype;
P.S         = S;
P.nvars     = nv;
P.nrows     = nr;
P.h0        = h0;
P.Ir        = Ir;
P.Jc        = Jc;
P.nnz       = numel(Ir);
P.groupCols = groupCols;
P.groupNz   = groupNz;
P.ncolour   = ncol;
P.density   = nnz(S)/numel(S);
end

% -------------------------------------------------------------------------
function J = do_jacobian(P, f, x)
x = x(:);
central = strcmp(P.type,'central');
if ~central
    v0 = f(x);
end

h    = P.h0 * max(1, abs(x));      % relative step, per variable
vals = zeros(P.nnz, 1);

for g = 1:P.ncolour
    cols = P.groupCols{g};
    if isempty(cols); continue; end
    k = P.groupNz{g};                      % nonzeros owned by this colour
    if isempty(k); continue; end

    xp = x; xp(cols) = xp(cols) + h(cols); % perturb the whole group at once
    if central
        xm = x; xm(cols) = xm(cols) - h(cols);
        dv = f(xp) - f(xm);
        den = 2*h(P.Jc(k));
    else
        dv = f(xp) - v0;
        den = h(P.Jc(k));
    end
    % each row hears from at most one column of this colour, so the finite
    % difference can be attributed unambiguously
    vals(k) = dv(P.Ir(k)) ./ den;
end

J = sparse(P.Ir, P.Jc, vals, P.nrows, P.nvars);
end

% -------------------------------------------------------------------------
function col = colour_columns(S)
% Greedy largest-first colouring of the column-intersection graph: two
% columns may share a colour only if they never both touch the same row.
nv  = size(S,2);
col = zeros(1,nv);
deg = full(sum(S,1));
[~,order] = sort(deg,'descend');
for q = order
    rws = find(S(:,q));
    if isempty(rws)
        col(q) = 1;                   % column affects nothing
        continue
    end
    used = col(any(S(rws,:),1));
    used = used(used > 0);
    c = 1;
    while any(used == c); c = c + 1; end
    col(q) = c;
end
end
