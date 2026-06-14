function [] = boundary(s,r)
% plot a circle, centre:(x,y) radius:r
theta1 = 0.944:0.01:pi-0.944;
theta2 = 0.944+pi:0.01:2*pi-0.944;
if s == 1
    Circle1=r*cos(theta);
    Circle2=r*sin(theta);
    c=[123,14,52];
    plot(Circle1,Circle2,'c','linewidth',1);
elseif s == 2
    Circle1=r*cos(theta1);
    Circle2=1/(2)*r*sin(theta1);
    c=[123,14,52];
    plot(Circle1,Circle2,'c','linewidth',2);
    hold on
    Circle3=r*cos(theta2);
    Circle4=1/(2)*r*sin(theta2);
    plot(Circle3,Circle4,'c','linewidth',2);
elseif s == 3
    ezplot('x.^2+(y-(x.^2).^(1/3)).^2-5');
axis equal
end
