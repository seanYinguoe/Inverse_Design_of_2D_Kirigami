function S = triangle_general_deployment(m, n, a, b, gamma, alpha)
%TRIANGLE_GENERAL_DEPLOYMENT  Rotating-units sheet of GENERAL triangular panels.
%
%   S = triangle_general_deployment(m, n, a, b, gamma, alpha)
%
%   The anisotropic generalisation of TRIANGLE_DEPLOYMENT: the rigid panels are
%   arbitrary triangles rather than equilateral ones, so the sheet expands by
%   two different principal stretches. a = b, gamma = 60 deg recovers rotating
%   equilateral triangles exactly.
%
%   INPUTS
%     m, n   : number of cells along the two lattice directions; one cell holds
%              one up-triangle and one down-triangle, so 2*m*n panels
%     a, b   : the two panel edge lengths enclosing gamma
%     gamma  : included angle [rad]
%     alpha  : PANEL rotation [rad], alpha = 0 compact. Half the hinge opening
%              angle theta of TRIANGLE_DEPLOYMENT, i.e. alpha = theta/2, so the
%              equilateral sheet is fully open at alpha = 60 deg.
%
%   OUTPUT  S, the struct shape PLOT_FORWARD consumes, plus the anisotropy
%   fields of ROTATING_UNIT_METRICS and the exact-overlap fields.
%
%   CONSTRUCTION
%     Let r1, r2, r3 be the panel corners relative to its centroid, counter-
%     clockwise. Up-panels turn by -alpha, down-panels (the point reflection,
%     corners -r_j) by +alpha. Requiring pinned corners to stay coincident,
%     the vector from an up-centroid to its three down-neighbours is
%
%        d_i(alpha) = R(-alpha) r_i + R(+alpha) r_{i-1},     i = 1,2,3
%
%     so the lattice is s = d1 - d2, t = d1 - d3, the down sublattice is the up
%     one offset by d1, and the pinning is
%
%        U(p,q) corner 1 <-> D(p,q)   corner 3
%        U(p,q) corner 2 <-> D(p-1,q) corner 1
%        U(p,q) corner 3 <-> D(p,q-1) corner 2
%
%     Each up-panel is pinned at all three corners to three different down-
%     panels and vice versa, so the mechanism closes for ANY triangle, not only
%     the equilateral one.
%
%   COLOURS
%     Six cuts meet at each rotation centre, so six colours are used - the same
%     scheme as TRIANGLE_DEPLOYMENT, each of the three edge families split by
%     the parity of its index along the lattice.
%
%   See also ROTATING_UNIT_METRICS, TRIANGLE_DEPLOYMENT, PLOT_FORWARD.

MET = rotating_unit_metrics('triangle', a, b, gamma, alpha);

% panel corners relative to the centroid, counter-clockwise
u = [a; 0];   v = b*[cos(gamma); sin(gamma)];
V = [[0;0], u, v];
if det([u v]) < 0, V = V(:,[1 3 2]); end
r = V - mean(V,2);

R  = @(ang) [cos(ang) -sin(ang); sin(ang) cos(ang)];
Rm = R(-alpha);   Rp = R(alpha);

d = zeros(2,3);
for i = 1:3, d(:,i) = Rm*r(:,i) + Rp*r(:,mod(i-2,3)+1); end
s = d(:,1) - d(:,2);
t = d(:,1) - d(:,3);

% ---- panels ---------------------------------------------------------------
Uv = Rm*r;                 % up-panel corners about its centroid
Dv = Rp*(-r);              % down-panel corners about its centroid
panels = cell(2*m*n,1);
idU = zeros(m,n);   idD = zeros(m,n);
k = 0;
for q = 0:m-1
    for p = 0:n-1
        c = p*s + q*t;
        k = k + 1;   panels{k} = (c + Uv).';          idU(q+1,p+1) = k;
        k = k + 1;   panels{k} = (c + d(:,1) + Dv).'; idD(q+1,p+1) = k;
    end
end
get = @(tbl,p,q) lookup(tbl, p, q, m, n);

% ---- shared (cut) edges and hinges ---------------------------------------
% {U edge, D offset [dp dq], D edge, colour, U pin local, D pin local}
pair = { [2 3], [ 0 -1], [2 3], @(p,q) 1 + mod(p,2), 3, 2
         [1 2], [-1  0], [1 2], @(p,q) 3 + mod(q,2), 2, 1
         [1 3], [ 0  0], [1 3], @(p,q) 5 + mod(q,2), 1, 3 };

seg = zeros(0,4);   cid = zeros(0,1);
hA  = zeros(0,2);   hB  = zeros(0,2);
for q = 0:m-1
    for p = 0:n-1
        iu = idU(q+1,p+1);
        for e = 1:size(pair,1)
            idd = get(idD, p+pair{e,2}(1), q+pair{e,2}(2));
            if idd == 0, continue; end
            eu = pair{e,1};   ed = pair{e,3};
            seg(end+1,:) = [panels{iu}(eu(1),:) panels{iu}(eu(2),:)];  %#ok<AGROW>
            cid(end+1,1) = pair{e,4}(p,q);                             %#ok<AGROW>
            seg(end+1,:) = [panels{idd}(ed(1),:) panels{idd}(ed(2),:)];%#ok<AGROW>
            cid(end+1,1) = pair{e,4}(p,q);                             %#ok<AGROW>
            hA(end+1,:)  = panels{iu }(pair{e,5},:);                   %#ok<AGROW>
            hB(end+1,:)  = panels{idd}(pair{e,6},:);                   %#ok<AGROW>
        end
    end
end

% ---- exact overlap test ---------------------------------------------------
PS = repmat(polyshape, numel(panels), 1);
for j = 1:numel(panels)
    PS(j) = polyshape(panels{j}, 'Simplify', true, 'KeepCollinearPoints', false);
end
PU = union(PS);
solid_area = numel(panels)*a*b*sin(gamma)/2;
union_area = area(PU);

% ---- centre on the origin -------------------------------------------------
allv = vertcat(panels{:});
ctr  = (min(allv,[],1) + max(allv,[],1))/2;
for j = 1:numel(panels), panels{j} = panels{j} - ctr; end
seg = seg - [ctr ctr];   hA = hA - ctr;   hB = hB - ctr;

% ---- pack -----------------------------------------------------------------
M45 = rotating_unit_metrics('triangle', a, b, gamma, MET.alpha_iso);

S.type       = 'triangle';
S.alpha      = alpha;
S.theta      = 2*alpha;
S.panels     = panels;
S.edges.seg  = seg;
S.edges.cid  = cid;
S.hinges     = (hA + hB)/2;
S.hinge_gap  = max([0; sqrt(sum((hA-hB).^2,2))]);
S.ncolors    = 6;
S.scale      = sqrt(a*b*sin(gamma)/2);
S.view_span  = max(n*norm(M45.s), m*norm(M45.t));
S.solid_area = solid_area;
S.union_area = union_area;
S.overlap    = solid_area - union_area;

S.mu = MET.mu;   S.dirs = MET.dirs;   S.lambda = MET.lambda;
S.nu = MET.nu;   S.p    = MET.p;      S.s = MET.s;   S.t = MET.t;
S.F  = MET.F;    S.cellangle = MET.cellangle;   S.k = MET.k;
end

% -------------------------------------------------------------------------
function k = lookup(tbl, p, q, m, n)
%LOOKUP  Panel id at lattice index (p,q), or 0 outside the patch.
k = 0;
if q < 0 || q > m-1 || p < 0 || p > n-1, return; end
k = tbl(q+1, p+1);
end
