function [Q ,AA, AV, OA, Kappa] = confmat_MG(C, Ct)
N = length(Ct);
Nc = max(Ct);

Q=zeros(Nc+1,Nc+1);
for i=1:Nc
    Q(1:Nc,i)=(histc(C(Ct==i),[1:Nc]))';
end
for i=1:Nc
    Q(Nc+1,i)=Q(i,i)/sum(Q(:,i));
    Q(i,Nc+1)=Q(i,i)/sum(Q(i,:));
end

AV = mean(Q(1:Nc,Nc+1));
AA = mean(Q(Nc+1,1:Nc));
OA=trace(Q)/sum(sum(Q(1:end-1,1:end-1)));
Q(Nc+2,1:Nc) = sum(Q(1:Nc,1:Nc),1);
Q(1:Nc,Nc+2) = sum(Q(1:Nc,1:Nc),2);
Kappa = (N*trace(Q) - Q(Nc+2,1:Nc)*Q(1:Nc,Nc+2))/(N^2 - Q(Nc+2,1:Nc)*Q(1:Nc,Nc+2));