function flag=test(~,bound,code)

flag=1;
[n,m]=size(code);

for i=1:n
    if code(i)<bound(i,1) || code(i)>bound(i,2)
        flag=0;
    end

%     y=funcable(code);
%     if y(1,1)<-295.6584 || y(1,2)<-233.3136 || y(1,3)<-233.3136 || y(1,4)<-295.6584
%         flag=0;
%     end

end
