function S = rotating_unit_deployment(shape, m, n, a, b, gamma, alpha)
%ROTATING_UNIT_DEPLOYMENT  Build a rotating-units sheet of either panel family.
%
%   S = rotating_unit_deployment(shape, m, n, a, b, gamma, alpha)
%
%     shape 'parallelogram' -> PARALLELOGRAM_DEPLOYMENT, four panels per vertex
%     shape 'triangle'      -> TRIANGLE_GENERAL_DEPLOYMENT, six panels per vertex
%
%   The two builders return the same struct shape, so PLOT_FORWARD and the
%   analysis code treat them identically. See ROTATING_UNIT_METRICS for the
%   continuum quantities and for what a, b, gamma and alpha mean.

switch lower(shape)
    case {'parallelogram','quad','rectangle','square'}
        S = parallelogram_deployment(m, n, a, b, gamma, alpha);
    case {'triangle','tri'}
        S = triangle_general_deployment(m, n, a, b, gamma, alpha);
    otherwise
        error('rotating_unit_deployment:shape', ...
            'shape must be ''parallelogram'' or ''triangle'', got ''%s''.', shape);
end
end
