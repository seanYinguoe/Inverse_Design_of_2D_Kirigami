function v = shape(code,XY,r)
%SHAPE Preset implicit equations. Vase uses signed r to select its arc.
x=XY(:,1);y=XY(:,2);
switch code
    case 1,v=hypot(x,y)-r;
    case 2,v=x.^2/r^2+y.^2/(r^2/4)-1;
    case 3,v=(r^2-x.^2)/3-(y-r).^2;
    case 4,v=y-target_profile(make_target('wavy'),x,'top');
    case 5,v=(x.^2+y.^2-r).^3-x.^2.*y.^3;
    otherwise,error('kirigami:Shape','Unknown target code.');
end
end
