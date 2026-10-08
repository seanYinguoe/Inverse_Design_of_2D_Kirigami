function M = parallelogram_metrics(a, b, gamma, alpha)
%PARALLELOGRAM_METRICS  Anisotropic kinematics of the rotating-parallelograms mechanism.
%
%   Thin wrapper kept for the existing call sites. The implementation now
%   lives in ROTATING_UNIT_METRICS, which covers triangular panels too; this
%   is exactly rotating_unit_metrics('parallelogram', ...).
%
%   M = parallelogram_metrics(a, b, gamma, alpha)
%
%   Closed-form continuum description of a rotating-units sheet whose rigid
%   panels are parallelograms of edge lengths a and b with included angle
%   gamma. Rotating squares is the special case a == b, gamma == pi/2, and it
%   is the ONLY isotropic member of the family.
%
%   INPUTS
%     a, b   : panel edge lengths
%     gamma  : included angle between the two panel edges [rad], 0 < gamma < pi
%     alpha  : PANEL rotation [rad], any size. Panels counter-rotate by
%              +/- alpha, so the hinge opening angle of the regular-case code
%              is theta = 2*alpha. alpha = 0 is compact.
%
%   OUTPUT  M, a struct. Fields that depend on alpha follow its size.
%     .mu        [mu1 mu2], the two anisotropy factors, mu1 >= 1 >= mu2,
%                mu1*mu2 = 1 exactly. Independent of alpha.
%     .H         (a^2+b^2)/(2*a*b*sin(gamma)) = (mu1+mu2)/2. H = 1 iff isotropic.
%     .dirs      2-by-2, principal directions as COLUMNS, dirs(:,i) belonging
%                to mu(i). Independent of alpha - the principal frame is fixed
%                along the whole deployment path.
%     .lambda    2-by-N principal stretches, lambda(i,:) = cos(alpha) + mu(i)*sin(alpha)
%     .nu        1-by-N Poisson's ratio, defined on the PRINCIPAL axes as the
%                logarithmic (true-strain) ratio
%                   nu = -dln(lambda2)/dln(lambda1)
%                i.e. minor direction responding to major. nu == -1 identically
%                for the isotropic case.
%     .p         1-by-N porosity, 1 - 1/(lambda1*lambda2). The solid area per
%                panel is fixed and the area per panel grows by det F, so this
%                generalises the regular-case p = 1 - lambda^(-2).
%     .s, .t     2-by-N deployed lattice vectors (columns per alpha)
%     .u, .v     2-by-1 reference panel edge vectors
%     .cellangle 1-by-N angle between s and t [rad]
%     .F         2-by-2-by-N effective deformation gradient, reference lattice
%                to deployed lattice
%
%   DERIVATION
%     Let u, v be the panel edge vectors and J the rotation by +pi/2. Requiring
%     that hinged corners stay coincident as alternate panels rotate by
%     +/- alpha gives the deployed lattice vectors
%
%        s(alpha) = u*cos(alpha) - (J*v)*sin(alpha)
%        t(alpha) = v*cos(alpha) + (J*u)*sin(alpha)
%
%     so with U = [u v] the deformation gradient F = [s t]*U^(-1) is
%
%        F(alpha) = cos(alpha)*I - sin(alpha)*M,        M = J*U*J*U^(-1).
%
%     M is symmetric - hence F is too, its principal frame is that of M and is
%     the SAME for every alpha. Writing the rows of U as r1, r2,
%
%        M = (1/det U) * [ -|r2|^2   r1.r2 ;  r1.r2   -|r1|^2 ]
%
%     whose determinant is the Gram determinant over (det U)^2, i.e. exactly 1,
%     and whose trace is -(a^2+b^2)/(a*b*sin(gamma)). Its eigenvalues are
%     therefore -mu and -1/mu with mu + 1/mu = (a^2+b^2)/(a*b*sin(gamma)), and
%     the principal stretches are lambda_i = cos(alpha) + mu_i*sin(alpha).
%
%   See also ROTATING_UNIT_METRICS, PARALLELOGRAM_DEPLOYMENT, FORWARD_ANISOTROPY.

if nargin < 4, alpha = 0; end
M = rotating_unit_metrics('parallelogram', a, b, gamma, alpha);
end
