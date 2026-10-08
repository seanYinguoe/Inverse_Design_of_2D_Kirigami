function [error_value,curve] = profile_error(boundary,target,degree)
%PROFILE_ERROR Area between fitted lower boundary and target, in design units^2.
% As in the original objective, align their vertical offsets at the left grip.
validateattributes(boundary,{'numeric'},{'2d','ncols',2,'real','finite'});
[x,~,group]=unique(boundary(:,1));y=accumarray(group,boundary(:,2),[],@mean);
assert(numel(x)>degree,'Too few distinct boundary points for the requested polynomial degree.');
w=target.clampx;
% The source objective excludes grip points and extrapolates the polynomial
% to both grips. Require samples on both sides; retain their extent in output.
assert(min(x)<0&&max(x)>0,'Boundary samples must span both sides of the sheet.');
[p,~,mu]=polyfit(x,y,degree);
shift=polyval(p,-w,[],mu)-target_profile(target,-w,'bottom');
f=@(q)polyval(p,q,[],mu)-shift;
error_value=integral(@(q)abs(f(q)-target_profile(target,q,'bottom')),-w,w, ...
    'AbsTol',1e-8,'RelTol',1e-6);
q=linspace(-w,w,201)';curve=[q f(q)];
end
