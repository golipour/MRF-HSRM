function IR=SRM_MG(ff , Q)
%% Read Picture --------------------
ff=mat2gray(ff)*256;
[N M V]=size(ff);
%% ------------------
R=1:M*N;
RR=reshape(R,N,M);   % initial regions
f=zeros(M*N,V);
f=double(reshape(ff,[],V));

S1(:,1)=reshape(RR(1:end-1,:),[],1);
S1(:,2)=S1(:,1)+1;
S2(:,1)=reshape(RR(:,1:end-1),1,[]);
S2(:,2)=S2(:,1)+N;
S=vertcat(S1,S2);
S(:,3)=max(abs(f(S(:,1),:)-f(S(:,2),:)),[],2);

[~,index]=sort(S(:,3));

%% ----------------------
ER=f;                  % initial regins mean
L=ones(length(f),1);          % region pixels number will be used for mean calculation
x=log(3*(M*N)^2);             %sigma prim
g = 256;
%% ----------------------
for i=1:length(S)
     a=S(index(i),1:2);
     
     r=R(a(1));
     %rr=a(1);
     while R(r) ~= r
         r = R(r);
         %R(rr)=r;
         %rr=r;
     end
     Rig(1)=r;
      
     r=R(a(2));
     %rr=a(2);
     while R(r) ~= r
         r = R(r);
         %R(rr)=r;
         %rr=r;
     end
     Rig(2)=r;
     
     r1=min(L(Rig(1)),g)*log(L(Rig(1))+1);
     r2=min(L(Rig(2)),g)*log(L(Rig(2))+1);
     b1=g*sqrt((x+r1)/(2*Q*L(Rig(1))));
     b2=g*sqrt((x+r2)/(2*Q*L(Rig(2))));
     P=abs(ER(Rig(1),:)-ER(Rig(2),:))-(b1+b2);
         
    if Rig(1) ~= Rig(2) && sum(P <= 0)==V
        [NR k]=min([R(a(1)) R(a(2))]);           % new Region
        NL=L(Rig(1))+L(Rig(2));                  % new region size
    
        ER(Rig(k),:)=ER(Rig(1),:)*(L(Rig(1))/NL)+ER(Rig(2),:)*(L(Rig(2))/NL);
        L(Rig(k))=NL;
        L(Rig(mod(k,2)+1))=0;
        R(a(mod(k,2)+1))=NR;
        R(Rig(mod(k,2)+1))=NR;
    end
end

while sum(R(R) == R) < M*N
    R=R(R);
end

RR=reshape(R,N,M);
%% Rigion ordering
i=1;
min2=0;
Max=max(RR(:));
while min2~=Max
    min1=min(RR(RR>min2));
    RR(RR==min1) = i;
    min2=min1;
    i=i+1;
end
IR=RR;
