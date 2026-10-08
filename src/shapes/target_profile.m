function y = target_profile(target,x,side)
%TARGET_PROFILE Top/bottom boundary shared by plots, constraints and COMSOL.
% Coordinates use the same design units as the kinematic pattern.
if nargin<3,side='bottom';end
assert(strcmp(target.type,'builtin'),'A mechanical target must be a built-in profile.');
r=target.r;
switch target.code
    case 1, upper=sqrt(r^2-x.^2);
    case 2, upper=0.5*sqrt(r^2-x.^2);
    case 3, upper=r-sqrt((r^2-x.^2)/3);
    case 4, upper=0.45*cos(0.8*pi*(x-1.25))+1.5;
    otherwise,error('kirigami:Profile','This target has no single-valued top/bottom profile.');
end
assert(isreal(upper),'Target profile queried outside its real-valued domain.');
if strcmp(side,'bottom'),y=-upper;elseif strcmp(side,'top'),y=upper;
else,error('kirigami:Side','Side must be top or bottom.');end
end
