close all
clear all
clc
% add all files
addpath(genpath('..'))




lambdas = [5000 1000 500 100 50 10 5 1];
noiselevels = [5:5:90];
psnrs = zeros(length(lambdas),length(noiselevels));

for iii = 1:length(lambdas),
    for jjj = 1:length(noiselevels),
        
        Verb   = 1;
        noiselevel = noiselevels(jjj);
        lambda = lambdas(iii);
        corrupted_image_name = sprintf('cameraman_pp_%d.png',noiselevel);
        B_Clean = double(sum(imread('cameraman_cle.png'),3));
        B_Corrupted =  double(sum(imread(corrupted_image_name),3));
        B_Clean = B_Clean/max(B_Clean(:));
        B_Corrupted = B_Corrupted/max(B_Corrupted(:));
        
        uexact = B_Clean;
        f = B_Corrupted;
        PSNR  =   snr(f, uexact);
        str = ['Input PSNR = ', num2str(PSNR)];
        disp(str)
        MaxIter=7;
        tic
        [y nois_ma2] = amf(f*255,39,0,0);
        PSNR  =   snr(y, uexact);
        str = ['PSNR = ', num2str(PSNR)];
        disp(str)
        
        E = y(:)/255 - f(:);
        [~, index] = sort(abs(E), 'descend');
        thres = max(abs(E(index(floor(1*noiselevel/100*512*512)))),1e-5);
        D = (abs(y/255-f)>thres);
        D = double(D);
        D1 = D;
        
        for i=1:MaxIter
            u = tvinpaint(f,lambda,D,[],[],[],[]);
            E = u(:) - f(:);
            [~, index] = sort(abs(E), 'descend');
            thres = abs(E(index(floor(1*noiselevel/100*512*512))));
            D = (abs(u-f)>thres);
            D = double(D);
            PSNR  =   snr(u, uexact);
            str = ['Iter = ',num2str(i), ', PSNR = ', num2str(PSNR)];
            disp(str)
        end
        toc
        
        PSNR  =   snr(u, uexact);
        str = ['Output PSNR = ', num2str(PSNR)];
        disp(str)
        psnrs(iii,jjj) = PSNR;
        psnrs
        
        %         figure;
        %         subplot(1,3,1); imshow(uexact,[]);title('Original','fontsize',13);
        %         subplot(1,3,2); imshow(f,[]); title('Corrupted','fontsize',13);
        %         subplot(1,3,3); imshow(u,[]); title('Recovered','fontsize',13);
        %         pause(2)
    end
end