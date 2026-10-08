function target = make_target(spec,varargin)
%MAKE_TARGET Named preset, implicit @(x,y) boundary, or N-by-2 closed polyline.
% Presets: circle, ellipse, vase, wavy, heart (old codes 1:5 also accepted).
% Size is the radius/shape parameter; ClampX is the grip half-width.
p = inputParser;
p.addParameter('Size',1,@(v) isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>0);
p.addParameter('ClampX',[],@(v) isempty(v)||(isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>0));
p.parse(varargin{:});
target = struct('type','','name','','code',[],'r',p.Results.Size, ...
    'clampx',p.Results.ClampX,'fun',[],'pts',[]);
names = {'circle','ellipse','vase','wavy','heart'};
if ischar(spec)||isstring(spec)
    spec = find(strcmp(validatestring(spec,names),names));
end
if isnumeric(spec)&&isscalar(spec)&&ismember(spec,1:5)
    target.type='builtin'; target.code=spec; target.name=names{spec};
    if ismember(spec,2:4)&&isempty(target.clampx),target.clampx=2.5;end
    if ismember(spec,[2 3])
        assert(target.r>=target.clampx,'Size must cover the grip half-width.');
    end
elseif isa(spec,'function_handle')
    target.type='implicit';target.name='custom';target.fun=spec;
elseif isnumeric(spec)&&size(spec,2)==2&&size(spec,1)>=3&&all(isfinite(spec(:)))
    target.type='curve';target.name='custom';
    if norm(spec(1,:)-spec(end,:))<1e-12,spec(end,:)=[];end
    assert(size(spec,1)>=3,'A closed curve needs at least three vertices.');
    target.pts=spec;
else
    error('kirigami:Target','Use a preset name, @(x,y), or a closed N-by-2 polyline.');
end
end
