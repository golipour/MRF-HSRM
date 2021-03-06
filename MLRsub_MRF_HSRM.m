%% MLRsub_MRF_HSRM by Meysam Golipour 2015, contact meysam.golipour@modares.ac.ir or miusay@gmail.com
% M. Golipour, H. Ghassemian, and F. Mirzapour, "Integrating Hierarchical Segmentation Maps
% With MRF Prior for Classification of Hyperspectral Images in a Bayesian Framework," IEEE
% Trans. on Geosci. Remote Sens., Agu.2015.

%%
clear;clc;
addpath '.\Lib\GCmex'
addpath '.\Lib\MLR'
addpath '.\Lib\Funcs'

% load the hyperspectral image
load('.\Data\Salinas')
HSI = salinas;

% load the ground truth 
load('.\Data\gtm')

% load the training samples
load('.\Data\Trains')
train = Trains(:,:,1);
%train = load('.\Data\train')
clear Trains

%% Parameters ---- (for more information about these parameters setting refere to the paper)
%PCA ---------------------------
Nf = 10;
%HSRM --------------------------
Q1 = exp(10);
Q2 = exp(6);
K = 10; % HSRM segmentation hierarchy levels
%Boarders Cost map -------------
alpha = 10; % alpha>>1 -> object-based classification even for an small beta value
%MRF ---------------------------
beta = 2;  % beta = 0 -> pixel-wise classifier

%% PCA -----------------------------
fprintf('PCA ... ');
[Nx Ny Nz] = size(HSI);
X = reshape(HSI,[],Nz);
X = PCA_MG(X,Nf);
Im = reshape(X,Nx,Ny,Nf);

fprintf('done!\n\n');

clear X 

%% HSRM Segmentation ------------------
fprintf('HSRM Segmentation ... ');

sz = size(Im);
f = int16(mat2gray(reshape(Im,[],Nf))*255);

tic
RR = FastHSRMCorr_MG_mex(f,sz,Q1,Q2,K);  % Save hyerarchical HSRM segmentation map in RR
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
fprintf('Boarder Extraction ... ');

CostMap = BoarderExNew_MG(RR,alpha);  % extracting fuzzy border/no-border map' from HSRM resualt
figure('Name', 'Fuzzy border/no-border map');imshow(CostMap);

fprintf('done!\n\n');

clear RR

%% MLRsub -------------------------------
fprintf('MLRsub Classification ... ');

X = reshape(HSI,[],Nz);

train_label = nonzeros(train);
train_data = X(train~=0,:);
[train_label,Ix] = sort(train_label);
train_data = train_data(Ix,:);

tic
[C , P] = MLR_MG(X,train_label, train_data,175);
elapsedTime = toc;

pp = reshape(P,Nx,Ny,[]);
cc = reshape(C,Nx,Ny);

Ct = gtm(:);
tr = train(:);
C = cc(:);
[Q ,AA, AV, OA, Kappa] = Confmat_MG(C(Ct~=0 & tr==0), Ct(Ct~=0 & tr==0));

figure('Name', 'MLRsub pix-wise classifier');imshow(label2rgb(double(cc)));

str = sprintf('don! elapsedTime = %.1fsec ', elapsedTime);fprintf(str);
str = sprintf('OA=%.2f \n\n',OA*100);fprintf(str);

clear HSI train_data

%% MLRsub-MLL  ------------------------------
fprintf('MLL Regularization ...\n');

tic
cc = MRF_MG(pp,beta,CostMap,'NE'); % MLRsub-MLL (non-Adaptive MRF)
elapsedTime = toc;

C = cc(:);
[Q ,AA, AV, OA, Kappa] = Confmat_MG(C(Ct~=0 & tr==0), Ct(Ct~=0 & tr==0));

figure('Name', 'MLRsub-MLL classifier');imshow(label2rgb(double(cc)));

str = sprintf('... don! elapsedTime = %.1fsec ',elapsedTime);fprintf(str);
str = sprintf('OA=%.2f \n\n',OA*100);fprintf(str);

%% MLRsub_MRF_HSRM --------------------------
fprintf('MRF-HSRM Regularization ...\n');

tic
cc = MRF_MG(pp,beta,CostMap,'E'); % MLRsub_MRF_HSRM (Adaptive MRF)
elapsedTime = toc;

C = cc(:);
[Q ,AA, AV, OA, Kappa] = Confmat_MG(C(Ct~=0 & tr==0), Ct(Ct~=0 & tr==0));

figure('Name', 'MLRsub_MRF_HSRM classifier');imshow(label2rgb(double(cc)));

str = sprintf('... don! elapsedTime = %.1fsec ',elapsedTime);fprintf(str);
str = sprintf('OA=%.2f \n\n',OA*100);fprintf(str);
