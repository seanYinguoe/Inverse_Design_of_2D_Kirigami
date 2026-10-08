function outfile = create_svg_tessellation(tessellation, t, filename, w, fillet_r)
% CREATE_SVG_TESSELLATION  Export a kirigami squares tessellation to an SVG file.
%
%   create_svg_tessellation(tessellation)
%   create_svg_tessellation(tessellation, t)
%   outfile = create_svg_tessellation(tessellation, t, filename)
%   outfile = create_svg_tessellation(tessellation, t, filename, w)
%   outfile = create_svg_tessellation(tessellation, t, filename, w, fillet_r)
%
%   The analytical kinematic model treats every cut as a zero-gap slit, so
%   neighbouring rigid squares meet exactly at the hinge points.  In a real
%   sheet the hinges are small ligaments: the slit does not reach the shared
%   corner but stops a short distance "t" away from it.  This function
%   reproduces that behaviour using the SAME gap construction as the COMSOL
%   LiveLink model (see livelink_rigid.m, lines 87-89).
%
%   Each rigid quadrant square pivots about ONE corner of every internal edge
%   (the hinge) while the opposite corner opens into the rotating-squares
%   hole.  Therefore a cut is shortened by "t" ONLY at its hinge (ligament)
%   end and runs fully to the open end - this gives the pinwheel cut pattern
%   rather than a uniform grid.  The hinge end of every edge is found
%   automatically by deploying the pattern a tiny amount with
%   tessellation_deployment and checking which shared corners stay joined.
%
%   The cut itself has zero thickness, so it is drawn as a plain line.
%
%   INPUTS
%     tessellation : m-by-n cell array, each cell a 16-by-2 list of points
%                    describing the 4 rigid quadrant squares of a unit
%                    (see create_unit.m).  Must be the undeployed pattern,
%                    e.g. tessellation_deployment(m,n,length,0).
%     t            : cut gap / ligament length (default 0.05, as in main.m).
%     filename     : output file name (default 'tessellation.svg').
%     w            : cut width (default 0). Each undeployed cut segment is
%                    treated as the centreline of a rectangular void of
%                    this width; w = 0 reproduces the original
%                    zero-thickness line. Collinear cut segments that meet
%                    end-to-end (no real ligament gap) are merged into a
%                    single continuous cut before drawing.
%     fillet_r     : corner fillet radius applied to the open (free) ends
%                    of each cut void (default w/4), automatically clamped
%                    so it never exceeds half of the rectangle's width or
%                    length. Any cut end that meets the sheet's outer
%                    boundary is left square (no fillet), since it
%                    terminates flush against the material edge rather
%                    than opening into free space.
%
%   OUTPUT
%     outfile      : full path of the written svg.
%
%   The svg is written to
%     <...>/kinematic tessellation optimization/output/

if nargin < 2 || isempty(t)
    t = 0.05;
end
if nargin < 3 || isempty(filename)
    filename = 'tessellation.svg';
end
if nargin < 4 || isempty(w)
    w = 0;
end
if nargin < 5 || isempty(fillet_r)
    fillet_r = w/4;
end

[m, n] = size(tessellation);
tolc = 1e-6;    % geometric coincidence tolerance

% index sets of the 4 rigid quadrant squares inside one unit (create_unit.m)
quad = {[1 2 3 4], [5 6 7 8], [9 10 11 12], [13 14 15 16]};

% --- a small deployment to discover the hinge (ligament) corners ---------
% the side length of a unit (difference of two known undeployed corners)
L = abs(tessellation{1,1}(3,1) - tessellation{1,1}(9,1));
deformed = tessellation_deployment(m, n, L, 0.1);  % tiny opening of the pattern

% --- collect every quadrant-square edge ----------------------------------
% for each edge keep the undeployed endpoints (for drawing) and the deployed
% endpoints (to test which corner stays joined).
allpts = [];
keys   = {};
data   = {};   % each entry: [A B dA dB] with A,B undeployed, dA,dB deployed
for i = 1:m
    for j = 1:n
        P  = tessellation{i,j};
        Pd = deformed{i,j};
        allpts = [allpts; P]; %#ok<AGROW>
        for q = 1:4
            idx = quad{q};
            for e = 1:4
                ia = idx(e);
                ib = idx(mod(e,4)+1);
                A  = P(ia,:);   B  = P(ib,:);
                dA = Pd(ia,:);  dB = Pd(ib,:);
                k = edge_key(A, B, tolc);
                pos = find(strcmp(keys, k), 1);
                if isempty(pos)
                    keys{end+1} = k;                       %#ok<AGROW>
                    data{end+1} = [A B dA dB];             %#ok<AGROW>
                else
                    data{pos} = [data{pos}; A B dA dB];    %#ok<AGROW>
                end
            end
        end
    end
end

% --- split into solid boundary edges and shortened internal cuts ---------
boundary = {};
cuts     = {};
for e = 1:numel(keys)
    rows = data{e};
    A = rows(1,1:2);
    B = rows(1,3:4);
    if size(rows,1) == 1
        % edge used by a single tile -> outer boundary of the sheet (solid)
        boundary{end+1} = [A B];                          %#ok<AGROW>
        continue;
    end
    % internal edge shared by two tiles: line up the deployed endpoints so
    % that both tiles refer to the same physical corner A and corner B.
    dA = zeros(size(rows,1),2);
    dB = zeros(size(rows,1),2);
    for r = 1:size(rows,1)
        if norm(rows(r,1:2) - A) < tolc
            dA(r,:) = rows(r,5:6);  dB(r,:) = rows(r,7:8);
        else
            dA(r,:) = rows(r,7:8);  dB(r,:) = rows(r,5:6);
        end
    end
    ligA = norm(dA(1,:) - dA(2,:)) < tolc;   % corner A stays joined -> hinge
    ligB = norm(dB(1,:) - dB(2,:)) < tolc;   % corner B stays joined -> hinge

    % shorten the cut by t only at the hinge (ligament) end(s)
    % (same construction as livelink_rigid.m)
    v  = B - A;
    al = t / norm(v);
    a2 = A;  b2 = B;
    if ligA, a2 = A + al*v; end
    if ligB, b2 = B - al*v; end
    if norm(a2 - b2) > tolc
        cuts{end+1} = [a2 b2];                            %#ok<AGROW>
    end
end

% neighbouring quadrant squares each contribute their own (already
% hinge-shortened) piece of what is physically one continuous slit; where
% two such pieces are collinear and meet with no real ligament gap, merge
% them into a single segment.
cuts = merge_collinear_cuts(cuts, tolc);

% --- model -> svg pixel mapping -----------------------------------------
xmin = min(allpts(:,1)); xmax = max(allpts(:,1));
ymin = min(allpts(:,2)); ymax = max(allpts(:,2));
scale  = 200;                          % pixels per unit length
margin = 20;                           % pixel border
W = (xmax - xmin)*scale + 2*margin;
H = (ymax - ymin)*scale + 2*margin;

% note: svg y points downwards, so flip y
map = @(p) [ (p(1) - xmin)*scale + margin, ...
             (ymax - p(2))*scale + margin ];

% --- determine output location ------------------------------------------
thisdir = fileparts(mfilename('fullpath'));         % .../function
outdir  = fullfile(fileparts(fileparts(thisdir)), 'results', 'exports');   % .../output
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
outfile = fullfile(outdir, filename);

% --- write the svg -------------------------------------------------------
fid = fopen(outfile, 'w');
if fid < 0
    error('create_svg_tessellation:cannotOpen', ...
          'Could not open %s for writing.', outfile);
end
fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(fid, ['<svg xmlns="http://www.w3.org/2000/svg" ', ...
              'width="%.2f" height="%.2f" viewBox="0 0 %.2f %.2f">\n'], ...
              W, H, W, H);

% sheet boundary (solid, no gap)
fprintf(fid, '  <g stroke="#000000" stroke-width="1.5" stroke-linecap="round">\n');
for e = 1:numel(boundary)
    p1 = map(boundary{e}(1:2));
    p2 = map(boundary{e}(3:4));
    fprintf(fid, '    <line x1="%.3f" y1="%.3f" x2="%.3f" y2="%.3f"/>\n', ...
            p1(1), p1(2), p2(1), p2(2));
end
fprintf(fid, '  </g>\n');

% cuts (shortened by t at the hinge end -> ligament), drawn either as
% zero-thickness lines (w == 0) or as void outlines of width w centred on
% the original cut segment (w > 0), filleted only at free ends - an end
% that lands on the sheet's outer boundary is left square.
if w <= 0
    fprintf(fid, '  <g stroke="#000000" stroke-width="1" stroke-linecap="round">\n');
    for e = 1:numel(cuts)
        p1 = map(cuts{e}(1:2));
        p2 = map(cuts{e}(3:4));
        fprintf(fid, '    <line x1="%.3f" y1="%.3f" x2="%.3f" y2="%.3f"/>\n', ...
                p1(1), p1(2), p2(1), p2(2));
    end
    fprintf(fid, '  </g>\n');
else
    wpx = w*scale;                          % cut width in pixels
    fpx = fillet_r*scale;                   % fillet radius in pixels
    fprintf(fid, '  <g fill="none" stroke="#000000" stroke-width="1">\n');
    for e = 1:numel(cuts)
        A = cuts{e}(1:2);  B = cuts{e}(3:4);
        p1  = map(A);
        p2  = map(B);
        mid = (p1 + p2)/2;
        Lpx = norm(p2 - p1);
        ang = atan2d(p2(2) - p1(2), p2(1) - p1(1));

        onBoundaryA = abs(A(1)-xmin)<tolc || abs(A(1)-xmax)<tolc || ...
                      abs(A(2)-ymin)<tolc || abs(A(2)-ymax)<tolc;
        onBoundaryB = abs(B(1)-xmin)<tolc || abs(B(1)-xmax)<tolc || ...
                      abs(B(2)-ymin)<tolc || abs(B(2)-ymax)<tolc;
        r0 = 0; if ~onBoundaryA, r0 = min([fpx, Lpx/2, wpx/2]); end
        r1 = 0; if ~onBoundaryB, r1 = min([fpx, Lpx/2, wpx/2]); end

        d = rounded_rect_path(Lpx/2, wpx/2, r0, r1);
        fprintf(fid, '    <path d="%s" transform="translate(%.3f %.3f) rotate(%.3f)"/>\n', ...
                d, mid(1), mid(2), ang);
    end
    fprintf(fid, '  </g>\n');
end

fprintf(fid, '</svg>\n');
fclose(fid);

fprintf('SVG written to: %s\n', outfile);
end

% -------------------------------------------------------------------------
function k = edge_key(a, b, tol)
% order-independent key for an edge so that a-b and b-a collapse to one edge
a = round(a/tol)*tol;
b = round(b/tol)*tol;
p = sortrows([a; b]);
k = sprintf('%.6f_%.6f__%.6f_%.6f', p(1,1), p(1,2), p(2,1), p(2,2));
end

% -------------------------------------------------------------------------
function out = merge_collinear_cuts(cuts, tol)
% merge axis-aligned cut segments that are collinear and touch end-to-end
% (gap < tol) into a single longer segment. Segments separated by a real
% ligament gap (on the order of t, always >> tol) are left untouched.
horizKeys = {}; horizIv = {};
vertKeys  = {}; vertIv  = {};
out = {};
for i = 1:numel(cuts)
    seg = cuts{i};
    A = seg(1:2); B = seg(3:4);
    if abs(A(2) - B(2)) < tol
        y = A(2);
        key = sprintf('%.6f', round(y/tol)*tol);
        pos = find(strcmp(horizKeys, key), 1);
        iv = sort([A(1) B(1)]);
        if isempty(pos)
            horizKeys{end+1} = key; horizIv{end+1} = iv;  %#ok<AGROW>
        else
            horizIv{pos} = [horizIv{pos}; iv];             %#ok<AGROW>
        end
    elseif abs(A(1) - B(1)) < tol
        x = A(1);
        key = sprintf('%.6f', round(x/tol)*tol);
        pos = find(strcmp(vertKeys, key), 1);
        iv = sort([A(2) B(2)]);
        if isempty(pos)
            vertKeys{end+1} = key; vertIv{end+1} = iv;     %#ok<AGROW>
        else
            vertIv{pos} = [vertIv{pos}; iv];                %#ok<AGROW>
        end
    else
        out{end+1} = seg;                                  %#ok<AGROW>
    end
end
for k = 1:numel(horizKeys)
    y = str2double(horizKeys{k});
    merged = merge_intervals(horizIv{k}, tol);
    for r = 1:size(merged,1)
        out{end+1} = [merged(r,1) y merged(r,2) y];        %#ok<AGROW>
    end
end
for k = 1:numel(vertKeys)
    x = str2double(vertKeys{k});
    merged = merge_intervals(vertIv{k}, tol);
    for r = 1:size(merged,1)
        out{end+1} = [x merged(r,1) x merged(r,2)];        %#ok<AGROW>
    end
end
end

% -------------------------------------------------------------------------
function merged = merge_intervals(iv, tol)
% collapse a set of [lo hi] intervals into their union, joining any pair
% whose gap is smaller than tol
iv = sortrows(iv, 1);
merged = iv(1,:);
for i = 2:size(iv,1)
    if iv(i,1) <= merged(end,2) + tol
        merged(end,2) = max(merged(end,2), iv(i,2));
    else
        merged(end+1,:) = iv(i,:);                          %#ok<AGROW>
    end
end
end

% -------------------------------------------------------------------------
function d = rounded_rect_path(hl, hw, r0, r1)
% path for a hl*2-by-hw*2 rectangle centred on the origin, with its long
% axis along x. r0 fillets the two corners at x = -hl (the A end), r1
% fillets the two corners at x = +hl (the B end); either may be 0 for a
% square end.
d = sprintf('M %.3f,%.3f ', -hl+r0, -hw);
d = [d, sprintf('L %.3f,%.3f ', hl-r1, -hw)];
if r1 > 0
    d = [d, sprintf('A %.3f,%.3f 0 0 1 %.3f,%.3f ', r1, r1, hl, -hw+r1)];
end
d = [d, sprintf('L %.3f,%.3f ', hl, hw-r1)];
if r1 > 0
    d = [d, sprintf('A %.3f,%.3f 0 0 1 %.3f,%.3f ', r1, r1, hl-r1, hw)];
end
d = [d, sprintf('L %.3f,%.3f ', -hl+r0, hw)];
if r0 > 0
    d = [d, sprintf('A %.3f,%.3f 0 0 1 %.3f,%.3f ', r0, r0, -hl, hw-r0)];
end
d = [d, sprintf('L %.3f,%.3f ', -hl, -hw+r0)];
if r0 > 0
    d = [d, sprintf('A %.3f,%.3f 0 0 1 %.3f,%.3f ', r0, r0, -hl+r0, -hw)];
end
d = [d, 'Z'];
end
