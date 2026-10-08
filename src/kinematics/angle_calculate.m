function theta = angle_calculate(nodes)
%ANGLE_CALCULATE  Unsigned interior angle at the MIDDLE node of three nodes.
%
%   theta = angle_calculate(nodes)
%
%   INPUT
%     nodes : 3-by-2, the (x,y) coordinates of three nodes [n1; n2; n3]
%
%   OUTPUT
%     theta : angle in radians between the segments n2->n1 and n2->n3,
%             i.e. the corner angle at n2. Always in [0, pi] - it carries no
%             sign, so the orientation of the corner has to be tested
%             separately with ifoverlapping.
%
%   Building block of the angle conditions Eqs. (5) and (10) of the paper
%   (sum of angles around an interior corner = 2*pi, resp. = pi).
%
%   NUMERICS  This uses the atan2 form
%
%       theta = atan2(|u x v|, u.v)
%
%   rather than the law of cosines with acos. The two agree to within
%   rounding for well-conditioned corners, but acos((a^2+b^2-c^2)/(2ab))
%   fails in two ways that matter inside an optimiser:
%     * rounding can push the argument just outside [-1,1], and MATLAB's
%       acos then returns a COMPLEX number, which propagates into fmincon's
%       constraint vector and derails the solve;
%     * its derivative is unbounded as the argument approaches +-1, i.e.
%       exactly at the flat and folded corners the optimiser has to pass
%       through, so finite-difference gradients there are meaningless.
%   The atan2 form is well conditioned over the whole range and returns 0
%   rather than NaN when two nodes coincide.

u = nodes(1,:) - nodes(2,:);   % n2 -> n1
v = nodes(3,:) - nodes(2,:);   % n2 -> n3

cross_z = u(1)*v(2) - u(2)*v(1);   % z-component of u x v, = |u||v|sin(theta)
dot_uv  = u(1)*v(1) + u(2)*v(2);   %                       = |u||v|cos(theta)

theta = atan2(abs(cross_z), dot_uv);
end
