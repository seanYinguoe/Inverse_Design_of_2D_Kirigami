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
%             i.e. the corner angle at n2, obtained from the law of cosines.
%             Always in [0, pi] - it carries no sign, so the orientation of
%             the corner has to be tested separately with ifoverlapping.
%
%   Building block of the angle conditions Eqs. (5) and (10) of the paper
%   (sum of angles around an interior corner = 2*pi, resp. = pi).

% Edges of the triangle n1-n2-n3
vector1 = nodes(2,:) - nodes(1,:);   % n1 -> n2
vector2 = nodes(3,:) - nodes(2,:);   % n2 -> n3
vector3 = nodes(1,:) - nodes(3,:);   % n3 -> n1

% Their lengths
mag1 = norm(vector1);
mag2 = norm(vector2);
mag3 = norm(vector3);

% Law of cosines on the side opposite to n2
theta = acos((mag1^2 + mag2^2 - mag3^2) / (2 * mag1 * mag2));
end
