function [] = plot_boundary(s,r)
% plot a circle, centre:(x,y) radius:r
% theta1 = 0.944:0.01:pi-0.944;
% theta2 = 0.944+pi:0.01:2*pi-0.944;
theta1 = 0:0.01:2*pi;
theta2 = acos(2.5/r):0.001:acos(-2.5/r);
theta3 = acos(2.5/r)+pi:0.001:acos(-2.5/r)+pi;
if s == 1  % circle
    Circle1=r*cos(theta1);
    Circle2=r*sin(theta1);
    plot(Circle1,Circle2,'color',[0.6350 0.0780 0.1840],'linewidth',1.5);
elseif s == 2 % ellipse
    Circle1=sqrt(1)*r*cos(theta3);
    Circle2=sqrt(1)*1/(sqrt(4))*r*sin(theta3);
    plot(Circle1,Circle2,'color',[0.6350 0.0780 0.1840],'linewidth',1.5);
    hold on
    Circle3=sqrt(1)*r*cos(theta2);
    Circle4=sqrt(1)*1/(sqrt(4))*r*sin(theta2);    
    plot(Circle3,Circle4,'color',[0.6350 0.0780 0.1840],'linewidth',1.5);
    plot([-2.5,-2.5],[sqrt(1)*1/(sqrt(4))*r*sin(acos(-2.5/r)),sqrt(1)*1/(sqrt(4))*r*sin(acos(2.5/r)+pi)],'color',[0.6350 0.0780 0.1840],'linewidth',1.5);
    plot([2.5,2.5],[sqrt(1)*1/(sqrt(4))*r*sin(acos(-2.5/r)+pi),sqrt(1)*1/(sqrt(4))*r*sin(acos(2.5/r))],'color',[0.6350 0.0780 0.1840],'linewidth',1.5);
elseif s == 3  % vase
    Circle1=sqrt(1)*r*cos(theta3);
    Circle2=sqrt(1)*1/(sqrt(3))*r*sin(theta3);
    plot(Circle1,Circle2+r,'color',[0.6350 0.0780 0.1840],'linewidth',1.5);
    hold on
    Circle3=sqrt(1)*r*cos(theta2);
    Circle4=sqrt(1)*1/(sqrt(3))*r*sin(theta2);    
    plot(Circle3,Circle4-r,'color',[0.6350 0.0780 0.1840],'linewidth',1.5);
    plot([-2.5,-2.5],[-1/(sqrt(3))*r*sin(acos(-2.5/r)+pi)-r,1/(sqrt(3))*r*sin(acos(2.5/r)+pi)+r],'color',[0.6350 0.0780 0.1840],'linewidth',1.5);
    plot([2.5,2.5],[-1/(sqrt(3))*r*sin(acos(-2.5/r)+pi)-r,1/(sqrt(3))*r*sin(acos(2.5/r)+pi)+r],'color',[0.6350 0.0780 0.1840],'linewidth',1.5);
elseif s == 4 % wavy
    x = [-2.5:0.1:2.5];
    y1 = 0.45*cos(0.8*pi*(x-(pi/(0.8*pi))))+1.5;
    y2 = -0.45*cos(0.8*pi*(x-(pi/(0.8*pi))))-1.5;
    plot(x,y1,'color',[0.6350 0.0780 0.1840],'linewidth',1.5);
    plot(x,y2,'color',[0.6350 0.0780 0.1840],'linewidth',1.5);
elseif s == 5 % heart;
    eqn = @(x,y) (x.^2 + y.^2 - r).^3 - x.^2 .* y.^3;
    fimplicit(eqn, 'color', [0.6350 0.0780 0.1840], 'linewidth', 1.5);
axis off
end
