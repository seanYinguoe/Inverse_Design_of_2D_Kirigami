function res = shape(s,nodes,r)
%SHAPE  Implicit equation of the target boundary shape.
%
%   res = shape(s,nodes,r)
%
%   Returns the value of the implicit function f(x,y) of the selected target
%   curve. res = 0 means the node lies exactly on the target boundary, so the
%   constraint files use res directly as an equality constraint (Eq. (7) of
%   the paper, "boundary shape condition"). For the closed convex shapes
%   res < 0 also means "inside", which is used as an inequality constraint.
%
%   INPUTS
%     s     : target shape selector
%               1 - circle       x^2 + y^2 = r^2
%               2 - ellipse      x^2/r^2 + y^2/(r/2)^2 = 1  (semi-axes r, r/2)
%               3 - vase         arc of a squashed circle offset by r; pass
%                                +r for the bottom arc and -r for the top one
%               4 - wavy         y = 0.3*cos(pi*x) + 1.2
%               5 - customised   heart, (x^2+y^2-r)^3 = x^2*y^3
%     nodes : [x y] coordinates of one node
%     r     : characteristic size of the shape (radius / semi-axis)
%
%   OUTPUT
%     res   : value of the implicit function at the node
%
%   NOTE  plot_boundary.m draws these same curves, but with its own hard-coded
%   parameters; if a curve is changed here it must be changed there as well.

if s == 1 % circle  x^2 + y^2 = r^2
    res = (nodes(1)^2 + nodes(2)^2)^0.5 - r;
elseif s == 2 % ellipse  x^2/r^2 + y^2/(1/4*r^2) = 1
    res = nodes(1)^2/(r^2) + nodes(2)^2/(1/4*r^2) - 1;
elseif s == 3 % vase: arc centred on (0,r); call with -r for the opposite arc
    %res = nodes(1)^2/3 + nodes(2)^2 - 2*r*nodes(2)+2/3*r^2;
    res = (1/3*r^2 - 1/3*nodes(1)^2) - (nodes(2)-r)^2;
elseif s == 4 % wavy cosine boundary
    res = nodes(2) - 0.3*cos(pi*nodes(1)) - 1.2;
elseif s == 5 % customized (heart curve)
    res = (nodes(1)^2+nodes(2)^2-r)^3 - nodes(1)^2*nodes(2)^3;
end
end
