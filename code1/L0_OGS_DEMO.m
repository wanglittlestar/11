clc;clear all;close all;

randn('seed',0);rand('seed',0);
 addpath('cimg','util') 
% Img = double(rgb2gray(imread('building_org.png')))/255;
Img = double(imread('pepper.png'))/255;
imshow(Img)
rank(Img)
if size(Img,3) > 1
    Img = rgb2gray(Img);
end
[row, col] = size(Img);

row = int2str(row);
col = int2str(col);
imageSize = [row 'x' col];
% R=7;
% [x,y] = meshgrid(-R:R,-R:R);
% K= double(x.^2 + y.^2 <= R^2);
% K = K/sum(K(:));
 K  =   fspecial('average',1); % For denoising
f1 = imfilter(Img,K,'circular');
% f1 = double(f1);
level = 0.9;
Bn =  imnoise(f1,'salt & pepper',level);
f = Bn;

O = ones(size(Img));
O(f == 1) = 0;
O(f == 0) = 0;
% Img = double(Img)/255; 
% f = f/255;
opts.lam       = [10];
opts.grpSz     =3; % OGS group size
opts.Nit       = 50000;
opts.Nit_inner = 5;
opts.tol       = 1/255;
opts.O        = O;

% main function
% lam=opts.lam;
% for j=1:length(lam)
%       out = L0_OGS_ADMM1(f,Img,K,opts);
%    end
%                    maxpsnr=max(max(peaksnr));
%                    [y1]=find(peaksnr==maxpsnr);
%                    lam(y1)
[out,peaksnr,k,t] = L0_OGS_ADMM1(f,Img,K,opts);

% figure;
% imwrite(out,"buildingOGSdeblurring10.png"),
% save buildingOGSdeblurring10 peaksnr t
% imwrite(out,"TVOGSdeblurring90.png")
% save TVOGSpepperdenosing90 peaksnr t