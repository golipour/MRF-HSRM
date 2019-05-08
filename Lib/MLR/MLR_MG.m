function [C, P] = MLR_MG(X,train_label,train_data,dind)
% dind <= Nb
%  number of classes
no_classes = max(train_label);
% size of image
% number of bands and number of samples
[l_all n_all] = size(X');
%%
[train_label,Ix] = sort(train_label);
train_inv = train_data(Ix,:);
%% ---------------------
train_label = train_label';
train_inv = train_inv';
n_train = length(train_label);
y = train_label;
[v,d]   = eig(train_inv*train_inv'/n_train);
d = diag(d);
%
dtrue = d(dind);
%%  subspace classifier
g = [];
g_all= [];

% number of classes

for k_iter = 1:no_classes
    index_k = train_label == k_iter;
    train_k = train_inv(:,index_k);
    n_k   = size(train_k,2);
    [v,d]   = eig(train_k*train_k'/n_k);
    d = diag(d);
    
    sub_sp = d<dtrue;
    sub(k_iter) = sum(sub_sp);
    tau = sub(k_iter);
    
    
    P = v(:,1:tau)*v(:,1:tau)'*train_inv;
    gall = zeros(1,n_train);
    
    for num_k = 1:n_train
        gall(num_k) = sqrt(P(:,num_k)'*P(:,num_k));
    end
    g = [g; gall];
    
    ggall=zeros(1,n_all);
    
    P_all = v(:,1:tau)*v(:,1:tau)'*X';
    
    for iter_all = 1:n_all
        ggall(iter_all) = sqrt(P_all(:,iter_all)'*P_all(:,iter_all));
    end
    g_all = [g_all;ggall];
end


% regularization parameter
lambda = eps;
w = subspace_classifier(-g,y,lambda);

% compute the probablity
P = subspace_mlogistic(w,-g_all);
[maxp,C] = max(P);
C = C';
P = P';

