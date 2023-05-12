clc;clear all;close all;

randn('seed',0);rand('seed',0);
 addpath('cimg','util') 

% I = double(imread('Cameraman.tif'))/255; 
% H = fspecial('average',9);
% Level = 0.1;
% B = imfilter(I,H,'circular','conv');
% Bn = imnoise(B,'salt & pepper',Level);
I = double(imread(sprintf('pepper.png')))/255;
% [u,s,v]=svds(I,10);
% I=u*s*v'
R=7;
[x,y] = meshgrid(-R:R,-R:R);
K= double(x.^2 + y.^2 <= R^2);
H = K/sum(K(:));
%  H = fspecial('average',1);
Level = 0.1;
%  B=conv2padded(I,H);
B = imfilter(I,H,'circular','conv');
imshow(B)
Bn = imnoise(B,'salt & pepper',Level);
psnr(Bn,I)
ssim(Bn,I)
%% Characteristic matrix
O = adaptive_median_detection(B,Bn);
% %% Initial evaluation index
IPSNR = 10*(log10(size(I,1)*size(I,2))-log10(norm(I(:) - Bn(:))^2));
ISNR = 20*log10(norm(I(:))/norm(I(:) - Bn(:)));
% ISSIM = ssim(I,Bn);
%% Initial valua
maxit = 5000;
espilon = 1/255;
gamma = 1.618;
mu = [200];
beta1= [3000];
beta2= [10];


%%



%% image denoise
  [U,SNR,PSNR,SSIM,i,t] = DetectionTVL1ADMM(I,Bn,H,O,maxit,espilon,mu,beta1,beta2,gamma);
 
% 
% figure(1);
% subplot(131); imshow(I,[]);
% subplot(132); imshow(Bn,[]);
% subplot(133); imshow(U,[]);
% 
% figure(2);
% subplot(131); imshow(I);
% subplot(132); imshow(Bn);
% subplot(133); imshow(U);


%% image denoise
% PSNR =zeros(length(mu),length(beta1),length(beta2));
% U = cell(length(mu),length(beta1),length(beta2));
% for i = 1:length(mu)
%     for j = 1:length(beta1)
%         for k = 1:length(beta2)
%             t = cputime;
%             [U{i,j,k},PSNR(i,j,k),SNR,SSIM,ii,t2] = DetectionTVL1ADMM(I,Bn,H,O,maxit,espilon,mu(i),beta1(j),beta2(k),gamma);
% %             t = cputime -t;
%         end
%     end
% end

% s = size(SNR);
% maxsnr = max(max(max(SNR)));
% Lax = find(SNR>=maxsnr);
% [o,p,q] = ind2sub(s,Lax);
% Lox_ax = [o,p,q];

% U = U{o,p,q};
% figure(1);
% subplot(131); imshow(I,[]);
% subplot(132); imshow(Bn,[]);
% subplot(133); imshow(U,[]);

% imwrite(U,"TVL1pepperdeblurring10.png")
%  save RankTVL1pepperdeblurring90 SNR SSIM t mu beta1 beta2 gamma U PSNR