function [] = plot_boundary(s,r)
%PLOT_BOUNDARY  Draw the target boundary shape on the current axes.
%
%   plot_boundary(s,r)
%
%   Companion of shape.m for visualisation: shape.m gives the implicit
%   equation the optimiser drives the boundary nodes onto, this function draws
%   the same curve in dark red so the deployed pattern can be compared with
%   its target.
%
%   INPUTS
%     s : target shape. Either a built-in code (1 circle | 2 ellipse |
%         3 vase | 4 wavy | 5 heart) or a target struct from make_target,
%         in which case any implicit curve or point list is drawn.
%     r : characteristic size of the shape (radius / semi-axis)
%
%   The clamped edges of the sheet are at x = +-2.5, i.e. half the width of
%   the default 5-by-5 unit grid. For shapes 2-4 only the top and bottom
%   boundaries are curved, so the arcs are drawn between those two x values
%   and closed with straight vertical lines.
%
%   NOTE  the curve parameters are hard-coded here and are NOT read from
%   shape.m; keep the two files in step when a target shape is changed.

target_color = [0.6350 0.0780 0.1840];   % dark red
half_width   = 2.5;                       % x position of the clamped edges

% --- general target from make_target -------------------------------------
if isstruct(s)
    hold on
    switch s.type
        case 'curve'
            P = [s.pts; s.pts(1,:)];
            plot(P(:,1), P(:,2), 'color', target_color, 'linewidth', 1.5);
        case 'implicit'
            h = fimplicit(@(x,y) s.fun(x,y), 'color', target_color, 'linewidth', 1.5);
            % keep the drawn range near the data rather than fimplicit's default
            ax = axis; h.XRange = ax(1:2); h.YRange = ax(3:4);
        case 'builtin'
            plot_boundary(s.code, s.r);
    end
    if ~isempty(s.clampx)
        yl = ylim;
        plot([-s.clampx -s.clampx], yl, 'color', target_color, 'linewidth', 1.5);
        plot([ s.clampx  s.clampx], yl, 'color', target_color, 'linewidth', 1.5);
    end
    return
end

theta_full  = 0:0.01:2*pi;                                  % full circle
theta_upper = acos(half_width/r):0.001:acos(-half_width/r); % top arc, x: +2.5 -> -2.5
theta_lower = theta_upper + pi;                             % bottom arc

if s == 1  % circle
    plot(r*cos(theta_full), r*sin(theta_full), ...
        'color',target_color,'linewidth',1.5);

elseif s == 2 % ellipse, semi-axes r and r/2
    plot(r*cos(theta_lower), 0.5*r*sin(theta_lower), ...
        'color',target_color,'linewidth',1.5);
    hold on
    plot(r*cos(theta_upper), 0.5*r*sin(theta_upper), ...
        'color',target_color,'linewidth',1.5);
    % straight clamped edges at x = -+2.5
    plot([-half_width,-half_width], ...
         [0.5*r*sin(acos(-half_width/r)), 0.5*r*sin(acos(half_width/r)+pi)], ...
        'color',target_color,'linewidth',1.5);
    plot([half_width,half_width], ...
         [0.5*r*sin(acos(-half_width/r)+pi), 0.5*r*sin(acos(half_width/r))], ...
        'color',target_color,'linewidth',1.5);

elseif s == 3  % vase: two arcs pushed apart by +-r (concave sides)
    plot(r*cos(theta_lower), 1/sqrt(3)*r*sin(theta_lower)+r, ...
        'color',target_color,'linewidth',1.5);
    hold on
    plot(r*cos(theta_upper), 1/sqrt(3)*r*sin(theta_upper)-r, ...
        'color',target_color,'linewidth',1.5);
    plot([-half_width,-half_width], ...
         [-1/sqrt(3)*r*sin(acos(-half_width/r)+pi)-r, 1/sqrt(3)*r*sin(acos(half_width/r)+pi)+r], ...
        'color',target_color,'linewidth',1.5);
    plot([half_width,half_width], ...
         [-1/sqrt(3)*r*sin(acos(-half_width/r)+pi)-r, 1/sqrt(3)*r*sin(acos(half_width/r)+pi)+r], ...
        'color',target_color,'linewidth',1.5);

elseif s == 4 % wavy cosine top and bottom edges
    x  = -half_width:0.1:half_width;
    y1 =  0.45*cos(0.8*pi*(x-(pi/(0.8*pi))))+1.5;
    y2 = -0.45*cos(0.8*pi*(x-(pi/(0.8*pi))))-1.5;
    plot(x,y1,'color',target_color,'linewidth',1.5);
    plot(x,y2,'color',target_color,'linewidth',1.5);

elseif s == 5 % heart
    eqn = @(x,y) (x.^2 + y.^2 - r).^3 - x.^2 .* y.^3;
    fimplicit(eqn, 'color', target_color, 'linewidth', 1.5);
    axis off
end
end
