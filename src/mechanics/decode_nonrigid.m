function nodes = decode_nonrigid(x,pattern)
%DECODE_NONRIGID Fill unique quarter-grid coordinates, then reflect both axes.
G=pattern.grid;G(pattern.free)=x;
q=cell(pattern.rows/2,pattern.cols/2);
for i=1:size(q,1)
    for j=1:size(q,2)
        ij=[2*i 2*j;2*i+1 2*j;2*i+1 2*j+1;2*i 2*j+1; ...
            2*i-1 2*j;2*i-1 2*j+1;2*i-1 2*j-1;2*i 2*j-1;2*i+1 2*j-1];
        for k=1:9,q{i,j}(k,:)=reshape(G(ij(k,1),ij(k,2),:),1,2);end
    end
end
nodes=mirror_nonrigid(q);
end
