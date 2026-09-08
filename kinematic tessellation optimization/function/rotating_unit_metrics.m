function M = rotating_unit_metrics(shape, a, b, gamma, alpha)
%ROTATING_UNIT_METRICS  Anisotropic kinematics of a rotating-units mechanism.
%
%   M = rotating_unit_metrics(shape, a, b, gamma, alpha)
%
%   One continuum description covering both panel families:
%
%     'parallelogram' - four panels meet at a vertex. Panel edges a, b with
%                       included angle gamma. Rotating SQUARES is a = b,
%                       gamma = 90 deg.
%     'triangle'      - six panels meet at a vertex. Panel is the triangle
%                       with two edges a, b enclosing gamma, so the third
%                       edge is c = sqrt(a^2 + b^2 - 2ab cos gamma). Rotating
%                       EQUILATERAL triangles is a = b, gamma = 60 deg.
%
%   INPUTS
%     shape : 'parallelogram' | 'triangle'
%     a, b  : the two panel edge lengths that enclose gamma
%     gamma : included angle [rad], 0 < gamma < pi
%     alpha : PANEL rotation [rad], any size. Alternate panels counter-rotate
%             by +/- alpha; alpha = 0 is compact. For the regular cases this is
%             HALF the hinge opening angle theta used by SQUARE_DEPLOYMENT and
%             TRIANGLE_DEPLOYMENT, i.e. alpha = theta/2.
%
%   OUTPUT  M, a struct; alpha-dependent fields follow the size of alpha.
%     .shape, .a, .b, .gamma, .c        the geometry as given (c only for triangles)
%     .mu        [mu1 mu2], mu1 >= mu2 > 0, the two anisotropy factors.
%                INDEPENDENT of alpha.
%     .k         mu1*mu2. 1 for parallelograms, 3 for triangles.
%     .H         (mu1+mu2)/2. Isotropic iff H = sqrt(k).
%     .dirs      2-by-2, principal directions as COLUMNS, dirs(:,i) belonging
%                to mu(i). INDEPENDENT of alpha - the principal frame is fixed
%                along the whole deployment path, for BOTH families.
%     .lambda    2-by-N principal stretches, cos(alpha) + mu(i)*sin(alpha)
%                Taken from trace(M) and det(M) analytically, so the isotropic
%                case comes out exactly degenerate.
%     .nu        1-by-N Poisson's ratio on the principal axes, defined as the
%                logarithmic ratio  nu = -dln(lambda2)/dln(lambda1).
%                SINGULAR at alpha = atan(mu1), where lambda1 is stationary and
%                the denominator vanishes: nu diverges there for an anisotropic
%                panel and is exactly -1 for an isotropic one. Values returned
%                near that rotation are meaningful only in the isotropic case.
%     .p         1-by-N porosity, 1 - 1/(lambda1*lambda2)
%     .s, .t     2-by-N deployed lattice vectors
%     .A0        2-by-2 reference (compact) lattice vectors as columns
%     .F         2-by-2-by-N deformation gradient, reference lattice to deployed
%     .M         2-by-2 the matrix below
%     .alpha_iso atan(sqrt(k)), the panel rotation at which an ISOTROPIC member
%                of this family is fully open: 45 deg for parallelograms,
%                60 deg for triangles. A convention for the end of the sweep,
%                NOT an admissibility limit - see the note below.
%
%   THE COMMON LAW
%     For both families the deformation gradient is affine in the rotation,
%
%        F(alpha) = cos(alpha)*I - sin(alpha)*M
%
%     with M SYMMETRIC - which is why F is symmetric and its principal frame
%     does not move as the sheet opens. M is built from the panel geometry as
%
%        parallelogram   M = J*U*J*U^(-1),           U = [u v]
%        triangle        M = 3*J*[r1 -r3]*A0^(-1),   r_i centroid-relative corners
%
%     and in BOTH cases
%
%        trace(M) = -(sum of squared panel edge lengths) / (2 * panel area)
%        det(M)   = k,   k = 1 (parallelogram), k = 3 (triangle)
%
%     so the eigenvalues of M are -mu1, -mu2 with mu1*mu2 = k and
%     mu1 + mu2 = -trace(M), and the principal stretches are
%
%        lambda_i(alpha) = cos(alpha) + mu_i*sin(alpha).
%
%     Isotropy needs mu1 = mu2 = sqrt(k), i.e. sum(l^2) = 2*sqrt(k)*(2*area).
%     For parallelograms that is a^2+b^2 >= 2ab*sin(gamma), equality only for
%     the square; for triangles it is a^2+b^2+c^2 >= 4*sqrt(3)*area, the
%     Weitzenboeck inequality, equality only for the equilateral triangle. Each
%     family therefore has exactly ONE isotropic member.
%
%   NOTE ON alpha_iso
%     lambda_i is maximised at alpha = atan(mu_i), so for an anisotropic panel
%     the two principal directions peak at DIFFERENT rotations and there is no
%     single natural "fully open" state. alpha_iso is the isotropic member's
%     value, used only as a sweep limit. The physical limit is set by panels
%     touching - test that with the .overlap field of the deployment functions.
%
%   See also PARALLELOGRAM_METRICS, ROTATING_UNIT_DEPLOYMENT,
%   PARALLELOGRAM_DEPLOYMENT, TRIANGLE_GENERAL_DEPLOYMENT.

if nargin < 5, alpha = 0; end
validateattributes(a,{'numeric'},{'scalar','positive','finite'},mfilename,'a');
validateattributes(b,{'numeric'},{'scalar','positive','finite'},mfilename,'b');
validateattributes(gamma,{'numeric'},{'scalar','>',0,'<',pi},mfilename,'gamma');

Jm = [0 -1; 1 0];                       % rotation by +pi/2
u  = [a; 0];
v  = b*[cos(gamma); sin(gamma)];
c  = NaN;

switch lower(shape)
    case {'parallelogram','quad','rectangle','square'}
        shape = 'parallelogram';
        A0 = [u v];                     % reference lattice = the panel edges
        Mm = Jm*A0*Jm/A0;
        edge2 = 2*(a^2 + b^2);          % sum of squared edges, four of them
        area  = a*b*sin(gamma);
        k     = 1;

    case {'triangle','tri'}
        shape = 'triangle';
        % panel corners, counter-clockwise, relative to the centroid
        V = [[0;0], u, v];
        if det([u v]) < 0, V = V(:,[1 3 2]); end
        r = V - mean(V,2);
        % compact lattice: A0 = [d1-d2, d1-d3] with d_i = r_i + r_{i-1}
        d0 = zeros(2,3);
        for i = 1:3, d0(:,i) = r(:,i) + r(:,mod(i-2,3)+1); end
        A0 = [d0(:,1)-d0(:,2), d0(:,1)-d0(:,3)];
        Mm = 3*Jm*[r(:,1) -r(:,3)]/A0;
        c     = sqrt(a^2 + b^2 - 2*a*b*cos(gamma));
        edge2 = a^2 + b^2 + c^2;
        area  = a*b*sin(gamma)/2;
        k     = 3;

    otherwise
        error('rotating_unit_metrics:shape', ...
            'shape must be ''parallelogram'' or ''triangle'', got ''%s''.', shape);
end

Mm = (Mm + Mm.')/2;                     % symmetric by construction; kill round-off

% Eigenvalues from the ANALYTIC trace and determinant rather than from eig.
% mu1*mu2 = k is exact, so this keeps the isotropic case (mu1 = mu2 = sqrt(k))
% bitwise degenerate. Taking them from eig instead leaves a 1-ulp gap, which
% the stationary point of lambda then amplifies into a meaningless nu.
Hc = -trace(Mm)/2;                      % = (mu1 + mu2)/2
disc = max(Hc^2 - k, 0);                % clamp: negative only by round-off
mu = [Hc + sqrt(disc), Hc - sqrt(disc)];

[V2, D2] = eig(Mm);                     % eig used for the DIRECTIONS only
[~, ix] = sort(-diag(D2).', 'descend');
V2 = V2(:,ix);
V2 = V2 ./ vecnorm(V2);
V2 = V2 .* sign(V2(1,:) + (V2(1,:)==0).*V2(2,:));   % stable sign

% ---- alpha-dependent quantities ------------------------------------------
al = alpha(:).';
ca = cos(al);   sa = sin(al);

lambda = ca + mu(:).*sa;
dl     = -sa + mu(:).*ca;
nu     = -(dl(2,:)./lambda(2,:)) ./ (dl(1,:)./lambda(1,:));

F = zeros(2,2,numel(al));
s = zeros(2,numel(al));   t = zeros(2,numel(al));
for j = 1:numel(al)
    F(:,:,j) = ca(j)*eye(2) - sa(j)*Mm;
    st = F(:,:,j)*A0;
    s(:,j) = st(:,1);   t(:,j) = st(:,2);
end

M = struct( ...
    'shape',     shape, ...
    'a', a, 'b', b, 'gamma', gamma, 'c', c, ...
    'mu',        mu, ...
    'k',         k, ...
    'H',         edge2/(2*(2*area)), ...
    'dirs',      V2, ...
    'lambda',    lambda, ...
    'nu',        nu, ...
    'p',         1 - 1./(lambda(1,:).*lambda(2,:)), ...
    's',         s, ...
    't',         t, ...
    'u',         A0(:,1), ...
    'v',         A0(:,2), ...
    'A0',        A0, ...
    'cellangle', atan2(s(1,:).*t(2,:) - s(2,:).*t(1,:), sum(s.*t,1)), ...
    'F',         F, ...
    'M',         Mm, ...
    'alpha_iso', atan(sqrt(k)));
end
