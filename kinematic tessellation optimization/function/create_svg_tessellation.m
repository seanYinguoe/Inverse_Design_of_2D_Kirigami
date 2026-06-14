function outfile = create_svg_tessellation(tessellation, t, filename)
% CREATE_SVG_TESSELLATION  Export a kirigami squares tessellation to an SVG file.
%
%   create_svg_tessellation(tessellation)
%   create_svg_tessellation(tessellation, t)
%   outfile = create_svg_tessellation(tessellation, t, filename)
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
outdir  = fullfile(fileparts(thisdir), 'output');   % .../output
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

% cuts (zero-thickness lines, shortened by t at the hinge end -> ligament)
fprintf(fid, '  <g stroke="#000000" stroke-width="1" stroke-linecap="round">\n');
for e = 1:numel(cuts)
    p1 = map(cuts{e}(1:2));
    p2 = map(cuts{e}(3:4));
    fprintf(fid, '    <line x1="%.3f" y1="%.3f" x2="%.3f" y2="%.3f"/>\n', ...
            p1(1), p1(2), p2(1), p2(2));
end
fprintf(fid, '  </g>\n');

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
