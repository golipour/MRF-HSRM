function Y = PCA_MG(X,Nf)
Nz = size(X,2);
X = 255*mat2gray(X);
Cx = cov(X);
[e,A]=eig(Cx);
Y = X*e(:,Nz-Nf+1:Nz);
end