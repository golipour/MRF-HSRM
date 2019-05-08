function Bmap = BoarderEx_MG(RR, alpha)
K = size(RR,3);
Gx = abs(imfilter(RR,[-1 0 1;-2 0 2;-1 0 1])) > 0;
Gy = abs(imfilter(RR,[-1 0 1;-2 0 2;-1 0 1]')) > 0;
Grad = mat2gray(Gx.^2+Gy.^2 > 0);
delta = abs(1-Grad);
g = (alpha-1)*delta(:,:,1) + sum(delta(:,:,2:K),3)/(K-1);
Bmap = g + 0;
end