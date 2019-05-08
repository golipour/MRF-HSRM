function train = trainsam(gtm,per,MinS,MaxS) 
[M,N] = size(gtm); 
g = gtm(:);
Max=max(g);
for i=1:Max
    s=sum(g==i);
    [~ , index]=sort(rand(s,1));
    a=zeros(s,1);
    c=round(per*s);
    if c<MinS 
        c=MinS;
    elseif c>MaxS
        c=MaxS;
    end
    b=zeros(s,1);b(1:c)=i;
    a(index)=b;
    g(g==i)=a;
end
train=reshape(g,M,N);
