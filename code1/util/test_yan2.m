close all
clear all
clc
% add all files  
addpath(genpath('..'))

ACWMF  = 0;
Verb   = 1;

noiselevel = 0;
noiselevel2 = 70;
lambda = 7000;

uexact = double(sum(imread('cameraman512.tif'),3));
uexact = uexact/max(uexact(:));

if noiselevel2 == 25
    if noiselevel == 5
        load Cam_5_25(14.41).mat
    elseif noiselevel == 10
        load Cam_10_25(14.5249).mat
    elseif noiselevel ==15
        load Cam_15_25(14.4735).mat
    end
elseif noiselevel2 == 40
    if noiselevel == 5
        load Cam_5_40(12.46).mat
    elseif noiselevel == 10
        load Cam_10_40(12.39).mat     
    elseif noiselevel == 15
        load Cam_15_40(12.43).mat
    end
elseif noiselevel2 == 30
    if noiselevel == 5
        load Cam_5_30(10.2841).mat
    elseif noiselevel == 10
        load Cam_10_30(10.2576).mat
    elseif noiselevel == 15
        load Cam_15_30(10.2145).mat
    end    
elseif noiselevel2 == 50
    if noiselevel == 5
        load Cam_5_50(8.0806).mat
    elseif noiselevel == 10
        load Cam_10_50(8.0933).mat
    elseif noiselevel == 15
        load Cam_15_50(8.0451).mat
    end
elseif noiselevel2 == 70
    if noiselevel == 5
        load Cam_5_70(6.6304).mat
    elseif noiselevel == 10
        load Cam_10_70(6.6026).mat
    elseif noiselevel == 15
        load Cam_15_70(6.5978).mat
    end
end

if noiselevel == 0 
    if noiselevel2 == 30
        load Cam_0_30.mat
    elseif noiselevel2 == 50
        load Cam_0_50.mat
    elseif noiselevel2 == 70
        load Cam_0_70.mat
    elseif noiselevel2 == 25
        load Cam_0_25.mat
    elseif noiselevel2 == 40
        load Cam_0_40.mat
    end
end

f3 = double(sum(imread('cameraman.png'),3));
f3 = f3/max(f3(:));
f      = f3;

%f      = imnoise(uexact,'salt & pepper',0.4);
%f      = addnoise(255*uexact,noiselevel2/100,'rd')/255;

if Verb
    PSNR  = snr(f,uexact);
    str = ['Input PSNR = ', num2str(PSNR)];
    disp(str)
end

MaxIter=7;
tic

if noiselevel2 == 25 || noiselevel2 == 40
    [y nois_ma2] = acwmf2(f*255);   
else
    [y nois_ma2] = amf(f*255,39,0,0);   
end

if Verb
    PSNR  = 20*log10(512/norm(y(:)/255-uexact(:)));
    str = ['ACWMF PSNR = ', num2str(PSNR)];
    disp(str)
end

E = y(:)/255 - f(:);
[~, index] = sort(abs(E), 'descend');
thres = max(abs(E(index(floor(1*noiselevel2/100*512*512)))),1e-5);
D = (abs(y/255-f)>thres);
D = double(D);
D1 = D;
for i=1:MaxIter
    u = tvinpaint(f,lambda,D,[],[],[],[]);
    E = u(:) - f(:);
    [~, index] = sort(abs(E), 'descend');
    thres = abs(E(index(floor(1*noiselevel2/100*512*512))));

    D = (abs(u-f)>thres);
    D = double(D);
    if noiselevel == 0
        lambda = lambda + 100;
    end
    if Verb
%         figure(5)
%         subplot(1,2,1)
%         imshow(D,[])
%         subplot(1,2,2)
%         imshow(u,[])
%         drawnow

        PSNR  = snr(u,uexact);
        str = ['Iter = ',num2str(i), ', PSNR = ', num2str(PSNR)];
        disp(str)
    end
end

if ACWMF
    [u, ~] = acwmf2(u*255);
    u      = double(u)/255;
end
toc

if Verb
    PSNR   = snr(u,uexact);
    str = ['Output PSNR = ', num2str(PSNR)];
    disp(str)
end

figure(1)
set(gcf,'Color',[1,1,1],'NumberTitle','off','Name','Exact vs Inpainted');
compareimages(uexact,'Exact',u,'Inpainted');
shg;
figure(2)
set(gcf,'Color',[1,1,1],'NumberTitle','off','Name','Input vs Inpainted');
compareimages(f,'Input',u,'Inpainted');
shg;
