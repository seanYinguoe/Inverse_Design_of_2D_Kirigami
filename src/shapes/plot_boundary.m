function plot_boundary(target,r)
%PLOT_BOUNDARY Draw the same boundary used by the optimiser.
if ~isstruct(target),target=make_target(target,'Size',r);end
hold on; color=[0.55 0.15 0.15];
if strcmp(target.type,'curve')
    P=[target.pts;target.pts(1,:)];plot(P(:,1),P(:,2),'Color',color);return
elseif strcmp(target.type,'implicit')
    fimplicit(target.fun,'Color',color);return
end
if target.code==1
    a=linspace(0,2*pi,400);plot(target.r*cos(a),target.r*sin(a),'Color',color);
elseif target.code==5
    fimplicit(@(x,y)(x.^2+y.^2-target.r).^3-x.^2.*y.^3,'Color',color);
else
    x=linspace(-target.clampx,target.clampx,300);
    bottom=target_profile(target,x,'bottom');top=target_profile(target,x,'top');
    plot([x fliplr(x) x(1)],[bottom fliplr(top) bottom(1)],'Color',color);
end
end
