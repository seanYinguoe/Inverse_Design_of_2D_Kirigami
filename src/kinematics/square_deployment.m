function S = square_deployment(m, n, l, theta)
%SQUARE_DEPLOYMENT  Forward kinematics of the rotating-SQUARES tessellation.
%
%   S = square_deployment(m, n, l, theta)
%
%   Thin wrapper around TESSELLATION_DEPLOYMENT that, besides the deformed
%   coordinates, also returns the connectivity needed to draw the mechanism:
%   which edges are cuts that open up, which colour each of them carries, and
%   where the hinges are.
%
%   INPUTS
%     m, n  : number of units along y and x (each unit = 2x2 quadrant squares,
%             so the mechanism is a 2m-by-2n checkerboard of rigid panels)
%     l     : side length of one unit (a quadrant square has side l/2)
%     theta : HINGE opening angle in radians, 0 = compact, pi/2 = fully open.
%             The panels themselves rotate by theta/2, which is the angle
%             TESSELLATION_DEPLOYMENT expects.
%
%   OUTPUT  S, see FORWARD_DEPLOYMENT for the field list.
%
%   PANEL AND EDGE BOOK-KEEPING
%     Quadrant square (I,J) of the 2m-by-2n checkerboard is stored with its
%     four corners in the local order [BL TL TR BR], so its local edges are
%
%        (1,2) left    (2,3) top    (3,4) right    (4,1) bottom
%
%     Two panels that touch along an edge in the compact state are pinned at
%     ONE end of that edge only; the edge is a cut everywhere else and opens
%     into a rhombic void as theta grows. Which end carries the pin
%     alternates like a checkerboard, mod(I+J,2), and reproduces exactly the
%     tie list documented in CREATE_UNIT.
%
%     Colours: a vertical cut is coloured by the parity of its row I, a
%     horizontal cut by the parity of its column J. The four cuts that meet
%     at any one rotation centre therefore always carry four different
%     colours, and a cut keeps its colour for every theta.

theta_max = pi/2;
if theta < -1e-9 || theta > theta_max + 1e-9
    warning('square_deployment:range', ...
        ['opening angle %.1f deg is outside the physical range [0, 90] deg ' ...
         'of the rotating-squares mechanism; panels will overlap.'], theta*180/pi);
end

tess = tessellation_deployment(m, n, l, theta/2);   % panels rotate by theta/2

nI = 2*m;   nJ = 2*n;      % checkerboard of rigid quadrant squares
s  = l/2;                  % side of one quadrant square

% ---- panels -------------------------------------------------------------
% rows of the 16-node unit list, and where each quadrant sits in the
% checkerboard:  I = 2*i + dI(q),  J = 2*j + dJ(q)
quad_rows = {1:4, 5:8, 9:12, 13:16};   % Q1 TR, Q2 BR, Q3 BL, Q4 TL
dI        = [ 0, -1, -1,  0];
dJ        = [ 0,  0, -1, -1];

panels = cell(nI*nJ,1);
id     = zeros(nI,nJ);
k = 0;
for i = 1:m
    for j = 1:n
        for q = 1:4
            k = k + 1;
            panels{k} = tess{i,j}(quad_rows{q},:);   % local order [BL TL TR BR]
            id(2*i+dI(q), 2*j+dJ(q)) = k;
        end
    end
end

% ---- shared (cut) edges and hinges --------------------------------------
nseg = 2*( nI*(nJ-1) + (nI-1)*nJ );
seg  = zeros(nseg,4);   cid = zeros(nseg,1);   e = 0;
nh   = nI*(nJ-1) + (nI-1)*nJ;
hA   = zeros(nh,2);     hB  = zeros(nh,2);     h = 0;

% vertical cuts: right edge of (I,J-1) against left edge of (I,J)
for I = 1:nI
    for J = 2:nJ
        A = panels{id(I,J-1)};   B = panels{id(I,J)};
        c = 3 + mod(I-1,2);
        seg(e+1,:) = [A(3,:) A(4,:)];   cid(e+1) = c;   % TR -> BR
        seg(e+2,:) = [B(1,:) B(2,:)];   cid(e+2) = c;   % BL -> TL
        e = e + 2;
        h = h + 1;
        if mod(I+J,2) == 0                              % pinned at the bottom
            hA(h,:) = A(4,:);   hB(h,:) = B(1,:);
        else                                            % pinned at the top
            hA(h,:) = A(3,:);   hB(h,:) = B(2,:);
        end
    end
end

% horizontal cuts: top edge of (I-1,J) against bottom edge of (I,J)
for I = 2:nI
    for J = 1:nJ
        A = panels{id(I-1,J)};   B = panels{id(I,J)};
        c = 1 + mod(J-1,2);
        seg(e+1,:) = [A(2,:) A(3,:)];   cid(e+1) = c;   % TL -> TR
        seg(e+2,:) = [B(4,:) B(1,:)];   cid(e+2) = c;   % BR -> BL
        e = e + 2;
        h = h + 1;
        if mod(I+J,2) == 0                              % pinned on the right
            hA(h,:) = A(3,:);   hB(h,:) = B(4,:);
        else                                            % pinned on the left
            hA(h,:) = A(2,:);   hB(h,:) = B(1,:);
        end
    end
end

% ---- pack ---------------------------------------------------------------
S.type      = 'square';
S.theta     = theta;
S.theta_max = theta_max;
S.panels    = panels;
S.edges.seg = seg;
S.edges.cid = cid;
S.hinges    = (hA + hB)/2;
S.hinge_gap = max([0; sqrt(sum((hA-hB).^2,2))]);   % 0 if the pinning is right
S.ncolors   = 4;
S.scale     = s;
S.lambda    = cos(theta/2) + sin(theta/2);         % linear expansion ratio
S.view_span = sqrt(2)*max(m,n)*l;                  % size of the fully open patch
end
