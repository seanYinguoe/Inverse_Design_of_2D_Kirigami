function S = triangle_deployment(m, n, a, theta)
%TRIANGLE_DEPLOYMENT  Forward kinematics of the rotating-TRIANGLES tessellation.
%
%   S = triangle_deployment(m, n, a, theta)
%
%   The triangular counterpart of TESSELLATION_DEPLOYMENT: rigid equilateral
%   triangles pinned to each other at their corners, i.e. the classical
%   rotating-triangles (twisted-kagome) auxetic mechanism with a single
%   degree of freedom.
%
%   INPUTS
%     m, n  : number of rhombic cells along the two lattice directions. One
%             cell holds one up-triangle and one down-triangle, so the patch
%             has 2*m*n rigid panels. Rows are offset so the compact patch
%             comes out roughly rectangular.
%     a     : side length of one triangle
%     theta : HINGE opening angle in radians. 0 = compact (the plain
%             triangular tiling), 2*pi/3 = 120 deg = fully open (the kagome
%             lattice, every hole a regular hexagon). Each panel itself
%             rotates by theta/2, up-triangles clockwise and down-triangles
%             counter-clockwise.
%
%   OUTPUT  S, see FORWARD_DEPLOYMENT for the field list.
%
%   KINEMATICS
%     Write the compact lattice as V(p,q) = p*u1 + q*u2 with u1 = [a 0] and
%     u2 = [a/2 a*sqrt(3)/2], and put
%
%        U(p,q) = [V(p,q)   V(p+1,q)   V(p,q+1)  ]     (up-triangle)
%        D(p,q) = [V(p+1,q) V(p,q+1)   V(p+1,q+1)]     (down-triangle)
%
%     Imposing that the three corners each triangle shares with a neighbour
%     stay coincident gives a one-parameter family in which the centroids
%     merely dilate about the origin,
%
%        lambda(theta) = 2*sin(theta/2 + pi/6),      lambda: 1 -> 2,
%
%     while each panel spins by -/+ theta/2 about its own centroid. The area
%     of the sheet therefore grows by lambda^2, i.e. by 4 at full opening, so
%     the mechanism has Poisson's ratio -1.
%
%   HINGES  (each triangle is pinned at all three of its corners)
%        U(p,q) 3 <-> D(p,q)   2      at V(p,q+1)
%        U(p,q) 2 <-> D(p,q-1) 3      at V(p+1,q)
%        U(p,q) 1 <-> D(p-1,q) 1      at V(p,q)
%
%   SHARED (CUT) EDGES AND THEIR COLOURS
%     The compact lattice has three edge directions; each is coloured by the
%     parity of its index along that direction, giving 2*3 = 6 colours. The
%     six cuts that meet at any one lattice vertex - the six edges that bound
%     the hexagonal void once the sheet opens - therefore always carry six
%     different colours, and a cut keeps its colour for every theta.
%
%        E1(p,q) = V(p,q)  -V(p+1,q)   U(p,q) 1-2  D(p,q-1) 2-3   1+mod(p,2)
%        E2(p,q) = V(p,q)  -V(p,q+1)   U(p,q) 1-3  D(p-1,q) 1-3   3+mod(q,2)
%        E3(p,q) = V(p+1,q)-V(p,q+1)   U(p,q) 2-3  D(p,q)   1-2   5+mod(q,2)

theta_max = 2*pi/3;
if theta < -1e-9 || theta > theta_max + 1e-9
    warning('triangle_deployment:range', ...
        ['opening angle %.1f deg is outside the physical range [0, 120] deg ' ...
         'of the rotating-triangles mechanism; panels will overlap.'], theta*180/pi);
end

u1  = [a 0];
u2  = [a/2 a*sqrt(3)/2];
lam = 2*sin(theta/2 + pi/6);          % dilation of the centroid lattice
Ru  = rot2(-theta/2);                 % up-triangles spin clockwise
Rd  = rot2( theta/2);                 % down-triangles counter-clockwise

% index window: row q uses p = p0(q) .. p0(q)+n-1. Successive rows are
% offset by u2, so pulling each one back by round(q/2) lattice steps stacks
% them like brickwork and keeps the compact patch rectangular instead of
% letting it shear away to the right.
p0   = @(q) -round(q/2);
pmin = p0(m-1);
pmax = n-1;
W    = pmax - pmin + 1;

idU = zeros(m,W);   idD = zeros(m,W);
panels = cell(2*m*n,1);
k = 0;
for q = 0:m-1
    for p = p0(q) + (0:n-1)
        V0  = p*u1 + q*u2;
        cU0 = V0 + [a/2, a/(2*sqrt(3))];        % centroid of U(p,q), compact
        cD0 = V0 + [a,   a/sqrt(3)    ];        % centroid of D(p,q), compact
        Uv0 = [V0; V0+u1; V0+u2];               % locals 1,2,3
        Dv0 = [V0+u1; V0+u2; V0+u1+u2];         % locals 1,2,3

        k = k + 1;
        panels{k} = lam*cU0 + (Uv0 - cU0)*Ru;
        idU(q+1, p-pmin+1) = k;

        k = k + 1;
        panels{k} = lam*cD0 + (Dv0 - cD0)*Rd;
        idD(q+1, p-pmin+1) = k;
    end
end

% ---- shared (cut) edges and hinges --------------------------------------
% pair table: {up local edge, down index offset [dp dq], down local edge,
%              colour, up hinge local, down hinge local}
pair = { [1 2], [ 0 -1], [2 3], @(p,q) 1+mod(p,2), 2, 3      % E1(p,q)
         [1 3], [-1  0], [1 3], @(p,q) 3+mod(q,2), 1, 1      % E2(p,q)
         [2 3], [ 0  0], [1 2], @(p,q) 5+mod(q,2), 3, 2 };   % E3(p,q)

seg = zeros(0,4);   cid = zeros(0,1);
hA  = zeros(0,2);   hB  = zeros(0,2);
for q = 0:m-1
    for p = p0(q) + (0:n-1)
        u = idU(q+1, p-pmin+1);
        for r = 1:size(pair,1)
            d = lookup(idD, p+pair{r,2}(1), q+pair{r,2}(2), pmin, m, W);
            if d == 0, continue; end          % neighbour outside the patch
            c  = pair{r,4}(p,q);
            eu = pair{r,1};   ed = pair{r,3};
            seg(end+1,:) = [panels{u}(eu(1),:) panels{u}(eu(2),:)];   %#ok<AGROW>
            cid(end+1,1) = c;                                         %#ok<AGROW>
            seg(end+1,:) = [panels{d}(ed(1),:) panels{d}(ed(2),:)];   %#ok<AGROW>
            cid(end+1,1) = c;                                         %#ok<AGROW>
            hA(end+1,:)  = panels{u}(pair{r,5},:);                    %#ok<AGROW>
            hB(end+1,:)  = panels{d}(pair{r,6},:);                    %#ok<AGROW>
        end
    end
end

% ---- pack ---------------------------------------------------------------
S.type      = 'triangle';
S.theta     = theta;
S.theta_max = theta_max;
S.panels    = panels;
S.edges.seg = seg;
S.edges.cid = cid;
S.hinges    = (hA + hB)/2;
S.hinge_gap = max([0; sqrt(sum((hA-hB).^2,2))]);
S.ncolors   = 6;
S.scale     = a;
S.lambda    = lam;
S.view_span = 2*max(n, m*sqrt(3)/2)*a;   % size of the fully open patch
end

% -------------------------------------------------------------------------
function R = rot2(t)
%ROT2  Rotation by t (counter-clockwise) acting on ROW vectors: v*rot2(t).
R = [cos(t) sin(t); -sin(t) cos(t)];
end

function k = lookup(tbl, p, q, pmin, m, W)
%LOOKUP  Panel id at lattice index (p,q), or 0 if it is outside the patch.
k = 0;
if q < 0 || q > m-1, return; end
c = p - pmin + 1;
if c < 1 || c > W,   return; end
k = tbl(q+1, c);
end
