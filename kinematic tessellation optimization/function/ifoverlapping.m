function res = ifoverlapping(nodes)
%IFOVERLAPPING  Signed orientation of a corner - the non-overlapping test.
%
%   res = ifoverlapping(nodes)
%
%   INPUT
%     nodes : 3-by-2, the (x,y) coordinates of three nodes [n1; n2; n3];
%             n2 is the corner, n1 and n3 the ends of the two adjacent edges.
%
%   OUTPUT
%     res : z-component of (n1-n2) x (n3-n2), i.e. <v1 x v2, n_hat> of
%           Eq. (6) of the paper with n_hat = (0,0,1).
%             res > 0 : the corner turns counter-clockwise (panels open)
%             res < 0 : the corner turns clockwise (panels have flipped over
%                       each other - overlap)
%             res = 0 : the three nodes are collinear
%           The magnitude is twice the triangle area, so it also measures how
%           far the corner is from collapsing.
%
%   The constraint files pass +res or -res to fmincon's inequality vector c
%   (fmincon enforces c <= 0), choosing the sign according to which way the
%   corner is expected to turn.

% Edges leaving the corner n2
vector1 = nodes(1,:) - nodes(2,:);
vector2 = nodes(3,:) - nodes(2,:);

% Cross product of the two edge vectors, lifted to 3D
cross_nodes = cross([vector1,0],[vector2,0]);

% Project on the out-of-plane unit vector n_hat = (0,0,1)
res = dot(cross_nodes,[0 0 1]);
end
