function cc = MRF_MG(pp,beta,Bmap,sw)
[Nx Ny Nc] = size(pp);
Dc = -log(pp+eps);
Sc = ones(Nc) - eye(Nc);     
% Expantion Algorithm ---------------
switch sw
    case 'E' 
        gch = GraphCut('open', Dc, beta*Sc, Bmap, Bmap);
    case 'NE' 
        gch = GraphCut('open', Dc, beta*Sc);
end
[gch seg] = GraphCut('expand',gch);
gch = GraphCut('close', gch);
% -----------------------------------
cc = seg+1;
end