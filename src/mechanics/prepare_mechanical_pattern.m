function pattern = prepare_mechanical_pattern(kinematic)
%PREPARE_MECHANICAL_PATTERN Check the hand-off; never silently mirror a different seed.
assert(kinematic.diagnostics.converged,'Kinematic constraints must converge first.');
T=kinematic.compact;[m,n]=size(T);mode=kinematic.config.mode;
assert(all(mod([m n],2)==0),'Mechanical symmetry needs even grid dimensions.');
assert(all(cellfun(@(p)isequal(size(p),[16 2]),T),'all'),'Expected 16 nodes per unit.');
P=vertcat(T{:});
assert(max(abs(imag(P)),[],'all')<1e-8&&all(isfinite(P),'all'),'Invalid compact geometry.');
P=real(P);offset=(min(P)+max(P))/2;
assert(max(abs((max(P)-min(P))-[n m]))<2e-3, ...
    'Use a fixed n-by-m sheet (FreeScale=false) for this mechanical parameterisation.');
T=cellfun(@(p)real(p)-offset,T,'UniformOutput',false);
pattern=struct('rows',m,'cols',n,'mode',mode,'target',kinematic.config.target);
assert(strcmp(pattern.target.type,'builtin')&&ismember(pattern.target.code,2:4), ...
    'Mechanical refinement currently supports ellipse, vase and wavy targets.');
assert(pattern.target.clampx>n/2,'Target grips must be outside the compact sheet.');
if strcmp(mode,'rigid')
    assert(m==2&&n==4,'The original rigid mechanical parameterisation is limited to 2-by-4.');
    ids=[1 2 3 4 5 6 8 9 10 14];
    q=cellfun(@(p)p(ids,:),T(1,1:2),'UniformOutput',false);
    P=vertcat(q{:});v=[P(:,1);P(:,2)];
    fraction=@(a,b,c)dot(a-b,c-b)/sum((c-b).^2);
    fractions=[fraction(P(4,:),P(7,:),P(3,:));fraction(P(19,:),P(7,:),P(3,:)); ...
        fraction(P(6,:),P(9,:),P(4,:));fraction(P(16,:),P(19,:),P(14,:))];
    pattern.x0=[v([5 2 15 12 7 3 29 34]);fractions]';
else
    ids=[1 2 3 4 5 8 9 13 14];
    G=zeros(m+1,n+1,2);counts=zeros(m+1,n+1);
    for i=1:m/2
        for j=1:n/2
            ij=[2*i 2*j;2*i+1 2*j;2*i+1 2*j+1;2*i 2*j+1; ...
                2*i-1 2*j;2*i-1 2*j+1;2*i-1 2*j-1;2*i 2*j-1;2*i+1 2*j-1];
            for k=1:9
                a=ij(k,1);b=ij(k,2);
                G(a,b,:)=G(a,b,:)+reshape(T{i,j}(ids(k),:),1,1,2);
                counts(a,b)=counts(a,b)+1;
            end
        end
    end
    G=G./counts;
    G(:,1,1)=-n/2;G(:,end,1)=0;G(1,:,2)=-m/2;G(end,:,2)=0;
    mask=true(size(G));mask(:,[1 end],1)=false;mask([1 end],:,2)=false;
    pattern.grid=G;pattern.free=find(mask);pattern.x0=G(mask)';
end
reconstructed=decode_pattern(pattern.x0,pattern);
expected=cellfun(@(p)p(ids,:),T,'UniformOutput',false);
pattern.reconstruction_error=max(abs(vertcat(reconstructed{:})-vertcat(expected{:})),[],'all');
assert(pattern.reconstruction_error<2e-3, ...
    'kirigami:UnrepresentableSeed', ...
    'Seed is not representable by the symmetric mechanical parameterisation (error %.3g).', ...
    pattern.reconstruction_error);
end
