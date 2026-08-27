function target = make_target(spec, varargin)
%MAKE_TARGET  Describe a target boundary shape, of any form.
%
%   target = make_target(code)                     one of the 5 built-ins
%   target = make_target(@(x,y) ..., 'Size', r)    an implicit curve
%   target = make_target(P)                        a closed curve, P is N-by-2
%   target = make_target(..., 'ClampX', w)         clamp the left/right edges
%
%   The program is not limited to the five hard-coded shapes. Eq. (7) of the
%   paper is already general - it asks that each boundary node p_i sit on the
%   target, measured by the distance to its projection onto that target - so
%   any closed curve works as long as we can evaluate "how far outside is
%   this point".
%
%   SPEC may be
%     a number 1..5  : the built-in shapes of shape.m
%                        1 circle  2 ellipse  3 vase  4 wavy  5 heart
%                      These keep their original hard-coded treatment so that
%                      published results still reproduce exactly.
%     a function handle @(x,y) -> value
%                      An implicit curve: the value must be 0 on the boundary,
%                      negative inside and positive outside. It has to be
%                      smooth, because the optimiser differentiates it.
%                      Example, a superellipse:
%                        f = @(x,y) (abs(x)/3).^4 + (abs(y)/2).^4 - 1;
%     an N-by-2 array  A closed curve through those points (the polygon is
%                      closed automatically; do not repeat the first point).
%                      The residual is the SIGNED DISTANCE to the polyline,
%                      negative inside - which is exactly the projection
%                      distance of Eq. (7). Use this for shapes you only have
%                      as data, e.g. traced from an image or exported from CAD.
%
%   NAME-VALUE OPTIONS
%     'Size'   : characteristic size r passed to the built-ins and available
%                to a function handle as target.r. Default 1.
%     'ClampX' : half-width w. When set, the left and right columns of nodes
%                are held straight at x = -+w instead of following the curve,
%                which is the tensile-grip boundary condition used for the
%                ellipse, vase and wavy cases. Default [] (follow the curve
%                all the way round).
%
%   OUTPUT
%     target : struct with fields type, code, fun, pts, r, clampx
%
%   EXAMPLES
%     t = make_target(1, 'Size', 3.05);                        % circle
%     t = make_target(@(x,y) x.^2/9 + y.^2/4 - 1, 'ClampX', 2.5);
%     th = linspace(0, 2*pi, 200).'; th(end) = [];
%     t = make_target([3*cos(th), 2*sin(th) + 0.4*sin(3*th)]); % free-form
%
%   Pass the result wherever the code takes a shape selector - rigid.m,
%   nonrigid.m, fit_initial_guess.m and design_freedom.m all accept it.

p = inputParser;
p.addParameter('Size', 1, @(v) isnumeric(v) && isscalar(v));
p.addParameter('ClampX', [], @(v) isempty(v) || (isnumeric(v) && isscalar(v)));
p.parse(varargin{:});

target = struct('type','', 'code',[], 'fun',[], 'pts',[], ...
                'r',p.Results.Size, 'clampx',p.Results.ClampX);

if isnumeric(spec) && isscalar(spec)
    if ~ismember(spec, 1:5)
        error('make_target:badCode', ...
              'Built-in shape code must be 1..5; got %g. Pass a function handle or a point list for anything else.', spec);
    end
    target.type = 'builtin';
    target.code = spec;

elseif isa(spec, 'function_handle')
    if nargin(spec) ~= 2
        error('make_target:badHandle', ...
              'An implicit curve must take two arguments, @(x,y).');
    end
    target.type = 'implicit';
    target.fun  = spec;

elseif isnumeric(spec) && ismatrix(spec) && size(spec,2) == 2 && size(spec,1) >= 3
    P = spec;
    if norm(P(1,:) - P(end,:)) < 1e-12
        P(end,:) = [];              % drop a repeated closing point
    end
    target.type = 'curve';
    target.pts  = P;

else
    error('make_target:badSpec', ...
        ['Target must be a code 1..5, a function handle @(x,y), or an ' ...
         'N-by-2 list of points defining a closed curve.']);
end
end
