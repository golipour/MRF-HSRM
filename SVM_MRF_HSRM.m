%% SVM_MRf_HSRM by Meysam Golipour 2015, contact meysam.golipour@modares.ac.ir or miusay@gmail.com
% M. Golipour, H. Ghassemian, and F. Mirzapour, "Integrating Hierarchical Segmentation Maps
% With MRF Prior for Classification of Hyperspectral Images in a Bayesian Framework," IEEE
% Trans. on Geosci. Remote Sens., Agu.2015.

%%
clear;clc;
addpath '.\Lib\libsvm-3.17';
addpath '.\Lib\libsvm-3.17\matlab';
addpath '.\Lib\libsvm-3.17\windows';
addpath '.\Lib\GCmex'
addpath '.\Lib\Funcs'

load('.\Data\gtm')
load('.\Data\HSI')
%load('.\Data\train') 
load('.\Data\Trains.mat')
train = Trains(:,:,1);

clear Trains

%% Parameters ---- (for more information about these parameters setting refere to the paper)
%PCA ---------------------------
Nf = 10;
%HSRM --------------------------
Q1 = exp(10);
Q2 = exp(6);
K = 10; % HSRM segmentation hierarchy levels
%Borders Cost map -------------
alpha = 10; % alpha>>1 -> object-based classification even for an small beta value
%MRF ---------------------------
beta = 10;  % beta = 0 -> pixel-wise classifier

%% PCA -----------------------------
fprintf('PCA ... ');
[Nx, Ny, Nz] = size(HSI);
X = reshape(HSI,[],Nz);
X = PCA_MG(X,Nf);
Im = reshape(X,Nx,Ny,Nf);
fprintf('done!\n\n');

clear X

%% HSRM Segmentation ------------------
fprintf('HSRM Segmentation ... ');

sz = size(HSI);
f = int16(mat2gray(reshape(Im,[],Nf))*255);

tic
RR = FastHSRMCorr_MG_mex(f,sz,Q1,Q2,K); % Save hyerarchical HSRM segmentation map in RR
elapsedTime = toc;

%{
% showing segmentation map hierarchy
for k = 1:K
 figure;imshow(label2rgb(RR(:,:,k), 'jet', 'W', 'shuffle'));   
end
%}

str = sprintf(' ... done! elapsedTime = %.1fsec\n\n',elapsedTime);fprintf(str);

clear Im f

%% Border Extraction ----------------
fprintf('Border Extraction ... ');

CostMap = BoarderExNew_MG(RR,alpha);  % extracting fuzzy border/no-border map' from HSRM resualt
figure('Name', 'Fuzzy border/no-border map');imshow(CostMap);

fprintf('done!\n\n');

clear RR

%% SVM -------------------------------
fprintf('SVM Classification ... \n');

% normalizing signature to optimize SVM result
X = reshape(HSI,[],Nz);
H = X; %preallocation of H for computational reasons
for i = 1:Nz
    m = min(X(train(:)~=0,i));
    M = max(X(train(:)~=0,i));
    H(:,i) = (X(:,i)-m)/(M-m);
end
X = 2*H-1;
%
train_label = nonzeros(train);
train_data = X(train~=0,:);
[train_label,Ix] = sort(train_label);
train_data = train_data(Ix,:);

tic
[C, P] = SVM_MG(X, train_data, train_label, gtm(:), 128, 0.0156);
elapsedTime = toc;

pp = reshape(P,Nx,Ny,[]);
cc = reshape(C,Nx,Ny);

Ct = gtm(:);
tr = train(:);
C = cc(:);
[Q ,AA, AV, OA, Kappa] = Confmat_MG(C(Ct~=0 & tr==0), Ct(Ct~=0 & tr==0));

figure('Name', 'SVM pix-wise classifier');imshow(label2rgb(double(cc)));

str = sprintf('... don! elapsedTime = %.1fsec ', elapsedTime);fprintf(str);
str = sprintf('OA=%.2f \n\n',OA*100);fprintf(str);

clear HSI train_data

%% SVM-MRF  ------------------------------
fprintf('MLL Regularization ... \n');

tic
cc = MRF_MG(pp,beta,CostMap,'NE'); % SVM-MLL(SVM-MRF) (non-Adaptive MRF)
elapsedTime = toc;

C = cc(:);
[Q ,AA, AV, OA, Kappa] = Confmat_MG(C(Ct~=0 & tr==0), Ct(Ct~=0 & tr==0));

figure('Name', 'SVM-MLL(SVM-MRF) classifier');imshow(label2rgb(double(cc)));

str = sprintf('... don! elapsedTime = %.1fsec ',elapsedTime);
fprintf(str);
str = sprintf('OA=%.2f \n\n',OA*100);
fprintf(str);

%% SVM-MRF-HSRM ----------------------------
fprintf('MRF-HSRM Regularization ... \n');

tic
cc = MRF_MG(pp,beta,CostMap,'E'); % SVM_MRF_HSRM (Adaptive MRF)
elapsedTime = toc;

C = cc(:);
[Q ,AA, AV, OA, Kappa] = Confmat_MG(C(Ct~=0 & tr==0), Ct(Ct~=0 & tr==0));

figure('Name', 'SVM_MRF_HSRM classifier');imshow(label2rgb(double(cc)));

str = sprintf('... don! elapsedTime = %.1fsec ',elapsedTime);
fprintf(str);
str = sprintf('OA=%.2f \n\n',OA*100);
fprintf(str);


