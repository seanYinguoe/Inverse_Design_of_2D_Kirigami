function res = shape(s,nodes,r)
% s = 1,circle; s = 2,ellispe; s = 3,heart.
if s == 1 % ellipse x^2/r^2 + y^2/(6*r^2) = 1  r = 3.2
    res = (nodes(1)^2 + nodes(2)^2)^0.5 - r;
elseif s == 2 % ellipse x^2/r^2 + y^2/(4*r^2) = 1  r = 3.2
    res = nodes(1)^2/(r^2) + nodes(2)^2/(1/4*r^2) - 1;
elseif s == 3 % vase r = 3.2
    %res = nodes(1)^2/3 + nodes(2)^2 - 2*r*nodes(2)+2/3*r^2;
    res = (1/3*r^2 - 1/3*nodes(1)^2) - (nodes(2)-r)^2;
elseif s == 4 % wavy
    res = nodes(2) - 0.3*cos(pi*nodes(1)) - 1.2;
elseif s == 5 % customized
    res = (nodes(1)^2+nodes(2)^2-r)^3 - nodes(1)^2*nodes(2)^3;
end
end
