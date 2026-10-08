function cmap = forward_colors()
%FORWARD_COLORS  Palette used to tag the shared (cut) edges of a tessellation.
%
%   cmap = forward_colors()
%
%   Returns a 6-by-3 RGB table. A square tessellation uses rows 1:4, a
%   triangular one uses all 6: that is exactly the number of shared edges
%   that meet at one rotation centre, so every hole of the deployed pattern
%   is surrounded by one edge of each colour.
%
%   Red is deliberately absent - it is reserved for the hinge markers.

cmap = [0.00 0.45 0.74;   % 1 blue
        0.93 0.53 0.00;   % 2 orange
        0.13 0.60 0.29;   % 3 green
        0.49 0.18 0.56;   % 4 purple
        0.00 0.68 0.75;   % 5 teal
        0.55 0.35 0.16];  % 6 brown
end
