  clc;clear all;close all;
randn('seed',0);rand('seed',0);
addpath('cimg','util')
opts.alpha=10;
opts.beta=0.1;   
opts.rho1=50;
opts.rho2=1;
miu=[34];
opts.iter=50000;         
lambdas= [0.1];
noiselevel=90; 
A=double(imread(sprintf('building.png')))/255;
[u,s,v]=svds(A,10);
B_Clean=(u*s*v')*255;
% B_Clean=double(imread(sprintf('building.png')));
P = GenBlurOper;
Amap = @(X)functionAX(P,X,'deblurring');
Atmap = @(X)functionAX(P',X,'deblurring');
B_Corrupted = impulsenoise(Amap(B_Clean),noiselevel/100,0);
imwrite(B_Corrupted/255,"lrbuildingdebluuring90.png") 
% Generate the mask matrix O
O = ones(size(B_Clean));
O(B_Corrupted==255)=0;
O(B_Corrupted==0)=0;
B_Clean = B_Clean/255;
B_Corrupted = B_Corrupted/255;
p =2;
P = GenBlurOper;
LargestEig = min(sqrt(sum(abs(P(:))>0)*sum(P(:).*P(:))), sum(abs(P(:))));% Largest Eigenvalue of A'A

acc =1/255;
pen_ratio =2;
 psnr=zeros(length(miu),length(lambdas));
        for j=1:length(lambdas)
                    for i=1:length(miu)
                            [U,SNR,psnr(i,j),iter,t2]= l0tv_padmm_color_for_plot_fobj11(B_Corrupted,O,Amap,Atmap,p,lambdas(j),LargestEig,acc,B_Clean,pen_ratio,opts,miu(i));
                         end
        end
                      maxpsnr=max(max(psnr));
                    [x1,y1]=find(psnr==maxpsnr);
                    miu(x1),lambdas(y1)
% [U,psnr(i,j)]= l0tv_padmm_color_for_plot_fobj11(B_Corrupted,O,Amap,Atmap,p,lambdas(j),LargestEig,acc,B_Clean,pen_ratio,opts,miu(i));
% 
%  psnr=zeros(length(lambdas));
%         for j=1:length(lambdas)
%              [U,SNR,psnr(j),iter,t2]=  l0tv_padmm_color_for_plot_fobj(B_Corrupted,O,Amap,Atmap,p,lambdas(j),LargestEig,acc,B_Clean,pen_ratio,opts);
%         end
%                     maxpsnr=max(max(psnr));
%                     [y1]=find(psnr==maxpsnr);
%                     lambdas(y1)
% [U,psnr,iter,t]= l0tv_padmm_color_for_plot_fobj11(B_Corrupted,O,Amap,Atmap,p,lambdas,LargestEig,acc,B_Clean,pen_ratio,opts,miu);
%   [U,psnr,iter,t]= l0tv_padmm_color_for_plot_fobj(B_Corrupted,O,Amap,Atmap,p,lambdas,LargestEig,acc,B_Clean,pen_ratio,opts);
% subplot(1,4,1); imshow(B_Clean,[]);title('Original','fontsize',13);
% subplot(1,4,2); imshow(B_Corrupted,[]); title('Corrupted','fontsize',13);
% subplot(1,4,3); imshow(U,[]); title('Recovered','fontsize',13);
% subplot(1,4,4); imshow(S,[]); title('Complement','fontsize',13);
% 
% imwrite(B_Corrupted,"housedeblurring30.png")
% % imwrite(U,"pepper10deblurring.png")
% 
% % subplot(1,4,3);
% imwrite(U," TVNL0housedeblurring90.png")
% save RankTVL0housedeblurring30 psnr t2 U SNR
% save RankTVNL0buildingdeblurring90 psnr t2 U SNR