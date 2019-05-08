function RR = HSRM_MG(ff, Q1, Q2, K)
ff = mat2gray(ff)*256;
[Nx, Ny, Nf] = size(ff);
R = 1:Nx*Ny;
RR = reshape(R,Nx,Ny);   % initial regions
f = zeros(Nx*Ny,Nf);
f = double(reshape(ff,[],Nf));

S1(:,1)=reshape(RR(1:end-1,:),[],1);
S1(:,2)=S1(:,1)+1;
S2(:,1)=reshape(RR(:,1:end-1),1,[]);
S2(:,2)=S2(:,1)+Nx;
S = vertcat(S1,S2);
S(:,3) = max(abs(f(S(:,1),:)-f(S(:,2),:)),[],2);
clear RR
[~,index] = sort(S(:,3));
% segments ----------------------
W2 = SRM_MG(ff,Q1);
clear ff
% ----------------------
R = reshape(W2,[],1);
MaxW = max(max(R)); 
ER = zeros(MaxW, Nf);         % initial regins mean
L = zeros(MaxW, 1);          % region pixels number will be used for mean calculation
for i = 1:MaxW
    %i
    ER(i, :) = mean(f(R==i, :));
    L(i) = sum(R==i);
end
% ----------------------
x = log(3*(Nx*Ny)^2);             %sigma prim
g = 256;
% ---------------------------------------
q = linspace(log(Q1),log(Q2),K);
iter = 1;
RR = zeros(Nx,Ny,length(q));
clear W2
for Q = exp(q) 
    Q
    Ss = S(R(S(:,1))~=R(S(:,2)),:);
    [~,index]=sort(Ss(:,3));
    for i = 1:length(Ss)
        a = Ss(index(i),1:2);
        Rig(1) = R(a(1));
        Rig(2) = R(a(2));
        if Rig(1) ~= Rig(2)
            r1 = min(L(Rig(1)),g)*log(L(Rig(1))+1);
            r2 = min(L(Rig(2)),g)*log(L(Rig(2))+1);
            b1 = g*sqrt((x+r1)/(2*Q*L(Rig(1))));
            b2 = g*sqrt((x+r2)/(2*Q*L(Rig(2))));
            P = abs(ER(Rig(1),:)-ER(Rig(2),:))-(b1+b2);
            
            if   sum(P <= 0)==Nf
                [NR k] = min(Rig);                         % new Region
                NL = L(Rig(1))+L(Rig(2));          % new region size
                
                ER(Rig(k),:) = ER(Rig(1),:)*(L(Rig(1))/NL)+ER(Rig(2),:)*(L(Rig(2))/NL);
                L(Rig(k)) = NL;
                L(Rig(mod(k,2)+1)) = 0;
                R(R==R(a(mod(k,2)+1))) = NR;    % update R
            end
        end
    end
    RR(:,:,iter) = reshape(R,Nx,Ny);
    iter = iter + 1;
end
end