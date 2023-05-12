clc;clear all;close all;

randn('seed',0);rand('seed',0);
 addpath('cimg','util') 
noiselevel=10;
opts.alpha =10;
opts.beta=0.1;  
opts.rho1=50; 
opts.rho2=1;
miu=200;
opts.iter=5000;
lambdas=[0.1];
B_Clean = double(imread(sprintf('building.png')));
B_Corrupted =impulsenoise(B_Clean,noiselevel/100,0)  ;
% Generate the mask matrix O
O = ones(size(B_Clean));
O(B_Corrupted==255)=0;
O(B_Corrupted==0)=0;
% 
B_Clean = B_Clean/255;
B_Corrupted = B_Corrupted/255;

p =2;
P = GenBlurOper;
LargestEig = min(sqrt(sum(abs(P(:))>0)*sum(P(:).*P(:))), sum(abs(P(:))));% Largest Eigenvalue of A'A

Amap = @(X)functionAX(P,X,'denoising');
Atmap = @(X)functionAX(P',X,'denoising');

acc = 1/255;pen_ratio = 10;
% % 
% t = cputime;
 psnr=zeros(length(miu),length(lambdas));
        for j=1:length(lambdas)
                    for i=1:length(miu)
                            [U,psnr(i,j),iter]= l0tv_padmm_color_for_plot_fobj11(B_Corrupted,O,Amap,Atmap,p,lambdas(j),LargestEig,acc,B_Clean,pen_ratio,opts,miu(i));
                        end
                    end
                    maxpsnr=max(max(psnr1));
                    [x1,y1]=find(psnr1==maxpsnr);
                    miu(x1),lambdas(y1)
% t=cputime-t;
%  [U,psnr(i,j)]= l0tv_padmm_color_for_plot_fobj11(B_Corrupted,O,Amap,Atmap,p,lambdas(j),LargestEig,acc,B_Clean,pen_ratio,opts,miu(i));
% %  
% psnr=zeros(length(lambdas));
%         for j=1:length(lambdas)
%              [U,psnr(j),iter,t]= l0tv_padmm_color_for_plot_fobj(B_Corrupted,O,Amap,Atmap,p,lambdas(j),LargestEig,acc,B_Clean,pen_ratio,opts);
%         end
%                     maxpsnr=max(max(psnr));
%                     [y1]=find(psnr==maxpsnr);
%                     lambdas(y1)
% % 
%  [U,psnr,iter,t]= l0tv_padmm_color_for_plot_fobj(B_Corrupted,O,Amap,Atmap,p,lambdas,LargestEig,acc,B_Clean,pen_ratio,opts);
%   [U,psnr,iter,t]= l0tv_padmm_color_for_plot_fobj11(B_Corrupted,O,Amap,Atmap,p,lambdas,LargestEig,acc,B_Clean,pen_ratio,opts,miu);

% S=ones(size(U)) - abs(B_Corrupted-U);
% subplot(1,4,1); imshow(B_Clean,[]);title('Original','fontsize',13);
% subplot(1,4,2); 
imshow(U,[]);
% subplot(1,4,3);
% imwrite(U," TVNL0pepperdenosing50.png")
% save TVL0pepperdenosing90 psnr t
% imshow(U,[]);;
% subplot(1,4,4); imshow(S,[]); title('Complement','fontsize',13);

