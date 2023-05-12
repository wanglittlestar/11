clc;clear all;close all;
randn('seed',0);rand('seed',0);

addpath(genpath('..'))

lambdas =  [0.2 0.4 0.6 0.8 1 2 4 6 8 10];
noiselevels = [10 20 30 40 50 60 70 80];
images = {'barbara','boat','cameraman','pirate','lenna'};
noisetypes = {'Salt-and-Pepper','Random-Valued'};
restorationtypes = {'denoising','deblurring'};

delete(sprintf('%s.txt',mfilename))

for mmm = 1:length(restorationtypes),
    for lll = 1:length(noisetypes),
        noisetype =noisetypes{lll};
        for kkk = 1:length(images),
            dlmwrite(sprintf('%s.txt',mfilename),datestr(now),'-append','delimiter','%s\t');
            mystr = sprintf('restorationtype:%s, noisetype:%s, image:%s',restorationtypes{mmm},noisetypes{lll},images{kkk});
            dlmwrite(sprintf('%s.txt',mfilename),mystr,'-append','delimiter','%s\t');
            psnrs = zeros(length(lambdas),length(noiselevels));
            psnrs2= zeros(length(lambdas),length(noiselevels));
            for iii = 1:length(lambdas),
                for jjj = 1:length(noiselevels),
                    lambda = lambdas(iii);
                    noiselevel = noiselevels(jjj);
                    B_Clean = double(sum(imread(sprintf('%s.png',images{kkk})),3));
                    
                    corrupted_image_name = sprintf('%s_%s_%s_%d.png',restorationtypes{mmm},images{kkk},noisetype,noiselevel);
                    B_Corrupted =  double(sum(imread(corrupted_image_name),3));
                    
                    % Generate the mask matrix O
                    O = ones(size(B_Clean));
                    if(strcmp(noisetype,'Salt-and-Pepper'))
                        max_val = max(abs(B_Corrupted(:))); min_val = min(abs(B_Corrupted(:))); 
                        O(B_Corrupted==max_val)=0; O(B_Corrupted==min_val)=0;
                    end
                    
                    B_Clean = B_Clean/max(B_Clean(:));
                    B_Corrupted = B_Corrupted/max(B_Corrupted(:));
                    
                    p = 2;
                    P = GenBlurOper;
                    LargestEig = min(sqrt(sum(abs(P(:))>0)*sum(P(:).*P(:))), sum(abs(P(:))));% Largest Eigenvalue of A'A
                    Amap = @(X)functionAX(P,X,restorationtypes{mmm});
                    Atmap = @(X)functionAX(P',X,restorationtypes{mmm});
                    
                    acc = 1/255;
                    [U] = l0tv_padmm_color(B_Corrupted,O,Amap,Atmap,p,lambda,LargestEig,acc,B_Clean);
                    [U2]  = l0tv_proj_reg(B_Corrupted,O,Amap,Atmap,p,lambda,LargestEig,acc,B_Clean);
                    
                    PSNR   =  snr(U, B_Clean);
                    PSNR2   =  snr(U2, B_Clean);
                    
                    psnrs(iii,jjj) = PSNR;
                    psnrs2(iii,jjj) = PSNR2;
                    psnrs
                    psnrs2
                    %     figure;
                    %     subplot(1,3,1); imshow(B_Clean,[]);title('Original','fontsize',13);
                    %     subplot(1,3,2); imshow(B_Corrupted,[]); title('Corrupted','fontsize',13);
                    %     subplot(1,3,3); imshow(U,[]); title('Recovered','fontsize',13);
                end
            end
            
            dlmwrite(sprintf('%s.txt',mfilename),'mpec-padmm','-append','delimiter','%s\t');
            dlmwrite(sprintf('%s.txt',mfilename),psnrs,'-append');
            dlmwrite(sprintf('%s.txt',mfilename),'optimal values','-append','delimiter','%s\t');
            dlmwrite(sprintf('%s.txt',mfilename),max(psnrs),'-append');

            dlmwrite(sprintf('%s.txt',mfilename),'penalty-decompositin-direct-projection','-append','delimiter','%s\t');
            dlmwrite(sprintf('%s.txt',mfilename),psnrs2,'-append');
            dlmwrite(sprintf('%s.txt',mfilename),'optimal values','-append','delimiter','%s\t');
            dlmwrite(sprintf('%s.txt',mfilename),max(psnrs2),'-append');
            
            dlmwrite (sprintf('%s.txt',mfilename),' ','-append','delimiter','%s\t')
            dlmwrite (sprintf('%s.txt',mfilename),' ','-append','delimiter','%s\t')
            dlmwrite (sprintf('%s.txt',mfilename),' ','-append','delimiter','%s\t')
        end
    end
end

