function res = shape_difference(p,s,r)
% s=1, ellipse; s=2, vase; s=3, wavy
if s == 1 % ellipse
    dis = -1/sqrt(4)*sqrt(1*r^2-(-2-0.5)^2) - polyval(p,-2-0.5);
    syms x1
    fun = abs((-1/sqrt(4)*sqrt(1*r^2-x1^2)-dis)-poly2sym(p,x1));
    Fint = int(fun,x1,-2.5,2.5);
    res = abs(Fint);
elseif s == 2 % vase
    dis = (sqrt(1/3*r^2-1/3*(-2-0.5)^2)-r) - polyval(p,-2-0.5);
    syms x1
    fun = abs((sqrt(1/3*r^2-1/3*(x1)^2)-r-dis) - poly2sym(p,x1));
    Fint = int(fun,x1,-2-0.5,2+0.5);
    res = vpa(Fint,3);
elseif s == 3 % wavy
    dis = -(0.2*cos(0.8*pi*(-2.5-1/0.8)) + 1.5) - polyval(p,-2-0.5);
    syms x1
    fun = abs(-(0.2*cos(0.8*pi*(x1-1/0.8)) + 1.5) - dis - poly2sym(p,x1));
    Fint = int(fun,x1,-2.5,2.5);
    res = vpa(Fint,3);
end
res = vpa(res,3);
end
