function v = target_value(target, XY)
%TARGET_VALUE  How far outside the target boundary is each point?
%
%   v = target_value(target, XY)
%
%   The one function every general target shape has to provide. Returns a
%   value per point that is
%       0   on the boundary
%     < 0   inside
%     > 0   outside
%   so the same number serves both jobs in the constraint files: as an
%   equality for the boundary nodes (Eq. 7 of the paper) and as an inequality
%   for every other node.
%
%   INPUTS
%     target : struct from make_target
%     XY     : N-by-2 array of query points
%
%   OUTPUT
%     v : N-by-1 values
%
%   For an implicit target this is just the curve's own function. For a target
%   supplied as a point list it is the SIGNED DISTANCE to the closed polyline,
%   which is precisely the projection distance ||p_i - p~_i|| of Eq. (7),
%   carrying a sign so that "inside" is usable as an inequality.
%
%   NOTE ON SMOOTHNESS  The signed distance to a polyline is continuous but
%   its derivative jumps where the nearest segment changes. fmincon's finite
%   differences tolerate that, but a densely sampled curve behaves better than
%   a coarse one, and a smooth implicit function better still. If you have a
%   formula for your shape, prefer the function-handle form.

switch target.type
    case 'implicit'
        v = target.fun(XY(:,1), XY(:,2));
        v = v(:);

    case 'curve'
        v = signed_distance(target.pts, XY);

    otherwise
        error('target_value:unsupported', ...
            ['target_value handles ''implicit'' and ''curve'' targets. ' ...
             'Built-in codes keep their own hard-coded treatment in ' ...
             'boundary_residual.m so that published results reproduce.']);
end
end

% -------------------------------------------------------------------------
function d = signed_distance(P, Q)
% Signed distance from each row of Q to the closed polygon P (negative inside).
A  = P;                 % segment starts
B  = P([2:end, 1], :);  % segment ends

nQ = size(Q,1);
d  = zeros(nQ,1);

AB   = B - A;                       % nP-by-2
LL   = sum(AB.^2, 2);               % squared segment lengths
LL(LL < eps) = eps;

for q = 1:nQ
    AQ = Q(q,:) - A;                            % nP-by-2
    t  = max(0, min(1, sum(AQ.*AB,2)./LL));     % projection parameter, clamped
    C  = A + t.*AB;                             % nearest point on each segment
    d(q) = sqrt(min(sum((Q(q,:) - C).^2, 2)));
end

in = inpolygon(Q(:,1), Q(:,2), P(:,1), P(:,2));
d(in) = -d(in);
end
