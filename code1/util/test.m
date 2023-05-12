function cpu_comparision

clc;clear all;close all;
randn('seed',0);rand('seed',0);
addpath('dataset');

load Cam_0_30.mat
n = numel(f3);

kkks = [round(n*0.7)];
% kkks = 7001;
for ikkk = 1:length(kkks),
    kkk = kkks(ikkk);
    type = 'denoising';
    B_Clean = double(sum(imread('cameraman512.tif'),3));
    B_Clean = B_Clean/max(B_Clean(:));
    
    p = 2;
    P = GenBlurOper;
    LargestEig = min(sqrt(sum(abs(P(:))>0)*sum(P(:).*P(:))), sum(abs(P(:))));% Largest Eigenvalue of A'A
    Amap = @(X)functionAX(P,X,type);
    Atmap = @(X)functionAX(P',X,type);
    
    B_Corrupted=f3;
    imshow(B_Corrupted)
    
    acc = 1e-8;
    [U] = l0tv2(B_Corrupted,Amap,Atmap,p,kkk,LargestEig,acc,B_Clean,P);
    imshow(U)
end




