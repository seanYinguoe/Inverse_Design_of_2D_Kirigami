function nodes = decode_pattern(x,pattern)
%DECODE_PATTERN Reconstruct the physical cut pattern from GA coordinates.
if strcmp(pattern.mode,'rigid')
    nodes=decode_rigid(x,pattern.rows,pattern.cols);
else
    nodes=decode_nonrigid(x,pattern);
end
assert(isreal(vertcat(nodes{:}))&&all(isfinite(vertcat(nodes{:})),'all'), ...
    'Candidate contains invalid coordinates.');
end
