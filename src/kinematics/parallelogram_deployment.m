function S = parallelogram_deployment(m, n, a, b, gamma, alpha)
%PARALLELOGRAM_DEPLOYMENT  Forward kinematics of a rotating-PARALLELOGRAMS sheet.
%
%   S = parallelogram_deployment(m, n, a, b, gamma, alpha)
%
%   The anisotropic generalisation of SQUARE_DEPLOYMENT: the rigid panels are
%   parallelograms of edge lengths a, b and included angle gamma instead of
%   squares, so the sheet expands by two different principal stretches. Setting
%   a == b and gamma == pi/2 recovers rotating squares exactly.
%
%   INPUTS
%     m, n   : number of panels along the v and u lattice directions
%     a, b   : panel edge lengths
%     gamma  : included angle between the panel edges [rad]
%     alpha  : PANEL rotation [rad]; alternate panels turn by +/- alpha.
%              alpha = 0 is the compact, as-cut state. Note this is HALF the
%              hinge opening angle theta used by SQUARE_DEPLOYMENT and
%              TRIANGLE_DEPLOYMENT: alpha = theta/2.
%
%   OUTPUT  S, the same struct shape SQUARE_DEPLOYMENT returns, so PLOT_FORWARD
%   draws it unchanged, plus the anisotropy fields:
%     .panels/.edges/.hinges/.hinge_gap/.ncolors/.scale/.view_span  as usual
%     .solid_area, .union_area
%                sum of the panel areas, and the area of their exact union.
%                They agree exactly while the panels do not overlap, so their
%                difference is an exact (not bounding-box) overlap test.
%
%   NOTE  No void polygons are returned. The panels of a rotating-units sheet
%   touch only at isolated hinge POINTS, so each void is pinched off from the
%   exterior by zero-width contacts and never registers as a hole of the panel
%   union - polyshape/HOLES returns nothing here. Extracting voids needs a
%   different construction (walk the four bounding panel edges), which is not
%   required by anything at present, so it is deliberately absent rather than
%   present and unverified.
%     .mu .dirs .lambda .nu .p .s .t .F .cellangle   from PARALLELOGRAM_METRICS
%
%   PANEL AND HINGE BOOK-KEEPING
%     Panel (I,J) of the m-by-n checkerboard sits at compact centre
%     (J-1/2)*u + (I-1/2)*v and turns by sigma*alpha with sigma = (-1)^(I+J).
%     Its corners are stored counter-clockwise as
%
%        1 = -(u+v)/2   2 = (u-v)/2   3 = (u+v)/2   4 = -(u-v)/2
%          (BL)            (BR)          (TR)          (TL)
%
%     so local edges are (1,2) bottom, (2,3) right, (3,4) top, (4,1) left.
%     Two panels sharing an edge are pinned at ONE end of it, and which end
%     alternates with sigma - that is what makes the sheet a mechanism rather
%     than a rigid tiling:
%
%        u-direction, A=(I,J) B=(I,J+1):  sigma_A>0 -> A2<->B1,  else A3<->B4
%        v-direction, A=(I,J) B=(I+1,J):  sigma_A>0 -> A3<->B2,  else A4<->B1
%
%     Colours follow SQUARE_DEPLOYMENT: a u-parallel cut is coloured by the
%     parity of its column J, a v-parallel cut by the parity of its row I, so
%     the four cuts meeting at any rotation centre always differ.
%
%   See also PARALLELOGRAM_METRICS, PLOT_FORWARD, SQUARE_DEPLOYMENT.

MET = parallelogram_metrics(a, b, gamma, alpha);

u = MET.u;   v = MET.v;                 % 2-by-1 reference edge vectors
s = MET.s;   t = MET.t;                 % 2-by-1 deployed lattice vectors

% corners of one panel, counter-clockwise, relative to its own centre
P0 = [-(u+v), (u-v), (u+v), -(u-v)]/2;  % 2-by-4, columns are corners 1..4

R = @(ang) [cos(ang) -sin(ang); sin(ang) cos(ang)];

% ---- panels ---------------------------------------------------------------
panels = cell(m*n,1);
id     = zeros(m,n);
k = 0;
for I = 1:m
    for J = 1:n
        sig = (-1)^(I+J);
        c   = (J-0.5)*s + (I-0.5)*t;    % deployed centre
        k = k + 1;
        panels{k} = (c + R(sig*alpha)*P0).';   % 4-by-2
        id(I,J)   = k;
    end
end

% ---- shared (cut) edges and hinges ---------------------------------------
nseg = 2*( m*(n-1) + (m-1)*n );
seg  = zeros(nseg,4);   cid = zeros(nseg,1);   e = 0;
nh   = m*(n-1) + (m-1)*n;
hA   = zeros(nh,2);     hB  = zeros(nh,2);     h = 0;

for I = 1:m                             % u-direction cuts: A right | B left
    for J = 2:n
        A = panels{id(I,J-1)};   B = panels{id(I,J)};
        c = 3 + mod(I-1,2);
        seg(e+1,:) = [A(2,:) A(3,:)];   cid(e+1) = c;
        seg(e+2,:) = [B(1,:) B(4,:)];   cid(e+2) = c;
        e = e + 2;   h = h + 1;
        if (-1)^(I+J-1) > 0             % sigma of the LEFT panel (I,J-1)
            hA(h,:) = A(2,:);   hB(h,:) = B(1,:);
        else
            hA(h,:) = A(3,:);   hB(h,:) = B(4,:);
        end
    end
end

for I = 2:m                             % v-direction cuts: A top | B bottom
    for J = 1:n
        A = panels{id(I-1,J)};   B = panels{id(I,J)};
        c = 1 + mod(J-1,2);
        seg(e+1,:) = [A(4,:) A(3,:)];   cid(e+1) = c;
        seg(e+2,:) = [B(1,:) B(2,:)];   cid(e+2) = c;
        e = e + 2;   h = h + 1;
        if (-1)^(I-1+J) > 0             % sigma of the LOWER panel (I-1,J)
            hA(h,:) = A(3,:);   hB(h,:) = B(2,:);
        else
            hA(h,:) = A(4,:);   hB(h,:) = B(1,:);
        end
    end
end

% ---- exact overlap test ---------------------------------------------------
% The panels of a mechanism touch only at hinge points, so while the motion is
% admissible the area of their union equals the sum of their areas exactly.
% Any shortfall is real overlap, measured by exact polygon boolean rather than
% by a bounding-box or corner-orientation proxy.
PS = repmat(polyshape, numel(panels), 1);
for k = 1:numel(panels)
    PS(k) = polyshape(panels{k}, 'Simplify', true, 'KeepCollinearPoints', false);
end
PU = union(PS);

solid_area = numel(panels)*a*b*sin(gamma);   % sum of the panel areas
union_area = area(PU);                       % exact, overlaps counted once

% ---- centre on the origin -------------------------------------------------
allv = vertcat(panels{:});
ctr  = (min(allv,[],1) + max(allv,[],1))/2;
for k = 1:numel(panels), panels{k} = panels{k} - ctr; end
seg = seg - [ctr ctr];   hA = hA - ctr;   hB = hB - ctr;

% ---- pack -----------------------------------------------------------------
s45 = MET.u*cos(pi/4) - [0 -1;1 0]*MET.v*sin(pi/4);
t45 = MET.v*cos(pi/4) + [0 -1;1 0]*MET.u*sin(pi/4);

S.type      = 'parallelogram';
S.alpha     = alpha;
S.theta     = 2*alpha;              % hinge angle, the regular-case convention
S.panels    = panels;
S.edges.seg = seg;
S.edges.cid = cid;
S.hinges    = (hA + hB)/2;
S.hinge_gap = max([0; sqrt(sum((hA-hB).^2,2))]);
S.solid_area = solid_area;
S.union_area = union_area;
S.overlap    = solid_area - union_area;   % > 0 exactly when panels intersect
S.ncolors   = 4;
S.scale     = sqrt(a*b);
S.view_span = max(n*norm(s45), m*norm(t45));

S.mu = MET.mu;   S.dirs = MET.dirs;   S.lambda = MET.lambda;
S.nu = MET.nu;   S.p    = MET.p;      S.s = MET.s;   S.t = MET.t;
S.F  = MET.F;    S.cellangle = MET.cellangle;
end
