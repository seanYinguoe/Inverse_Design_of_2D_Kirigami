clear;
clc;
% initial configuration of kirigami tessellation
m = 4;  % the number of tessellation units in each axis 
n = 4;
length = 1;  % the length of each units
Rotate_angels = pi/3; % the rotation angles of squares tessllation
tessellation_intial = tessellation_deployment(m,n,length,0);
tessellation_transformed = tessellation_deployment(m,n,length,Rotate_angels);
% % set the boudary condition
s = 1; % s = 1,circle; s = 2,ellispe; s = 3, vase; s = 4, wavy; s = 5, heart
r = 3.2*sqrt(2.0/2.2); % the radius of circle; ellipse: r = 3.2*sqrt(2.0/2.2);
p = 1; % p = 1, rigid; p = 2, nonrigid;
t = 0.30*1/m; % cut gap
w = 0.50*t; % cut width
% export the initial tessellation to an svg file (cuts shortened by t so the
% hinge points stay connected by a small ligament, drawn with width w and
% filleted corners)
create_svg_tessellation(tessellation_intial, t, 'tessellation_initial_8X8.svg', w);
% optimize the kirigami tessellation
%tessellation_optimized = tessellation_optimization(tessellation_transformed,s,r,p);
% compact the kirigami tessellation
%tessellation_compacted = tessellation_compaction(tessellation_optimized);

% plot result
figure(2);
subplot(2,2,1);                                  
title('intial kirigami tessellation');
plot_tessellation(tessellation_intial);
axis off

subplot(2,2,2);  
title('deployed kirigami tessellation');
plot_tessellation(tessellation_transformed);
plot_boundary(s,r);
axis off

subplot(2,2,3);
title('optimized kirigami tessellation');
plot_tessellation(tessellation_optimized);
% plot([-4:0.1:4],0.3*cos(pi*[-4:0.1:4]) + 1.3);
% plot([-4:0.1:4],-0.3*cos(pi*[-4:0.1:4]) - 1.3);
plot_boundary(s,r);
axis off

subplot(2,2,4);
title('valid compacted kirigami tessellation');
plot_tessellation(tessellation_compacted);
axis off