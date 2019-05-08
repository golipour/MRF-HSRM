function CC = MaxVot_MG(cc,RR)
R = RR(:);
C = cc(:);
Nc = max(C);
labels = unique(R);
CR = zeros(length(labels),1); 
i = 1;
for r = labels' 
    CR(i) = mode(C(R==r));
    i = i+1;
end
% parfor j = 1:length(labels); 
%     CR(j) = mode(C(R==labels(j)));
% end
CC = changem(RR,CR,labels);
end


