% TV inpainting demo using tvrestore, Pascal Getreuer 2010
%
% This is a MATLAB script.  To run it, enter on the MATLAB console
%   >> tvinpaint_demo

path(path,'C:\Users\yanm\Dropbox\BlindInpaiting\DECOLOR_MATLAB_WINDOWS\DECOLOR_MATLAB_WINDOWS\gco-v3.0\matlab')
path(path,'C:\Users\yanm\Dropbox\BlindInpaiting\DECOLOR_MATLAB_WINDOWS\DECOLOR_MATLAB_WINDOWS\internal')

%lambda = 1.25e1;
lambda = 1.25e1;
%lambda = 1.25;
uexact = double(sum(imread('cameraman.tif'),3))/255;
%uexact = phantom(512);
% D1 = (uexact == 1);
% uexact(D1) = 0;
uexact = uexact/max(uexact(:));
uexact1 = uexact;
uexact = uexact+10/255*randn(size(uexact));
uexact = imnoise(uexact,'salt & pepper', 0.1);

% Construct inpainting region
[x,y] = meshgrid(1:size(uexact,2),1:size(uexact,1));
[th,r] = cart2pol(x-size(uexact,2)/2,y-size(uexact,1)/2);
%D = (sin(r/1.5+th) > 0.950);
D = (sin(r/20.5+th) > 0.997);

f = uexact;
tic
%f = f3;
f(D) = 1*rand(nnz(D),1);
%f = f3;
% Make a figure 
clf;
set(gcf,'Color',[1,1,1],'NumberTitle','off','Name','TV Inpainting');
compareimages(f,'Input',f,'Inpainted');
shg;

% Inpaint
D = zeros(size(uexact,1),size(uexact,2));
D2 =D;
thres = 0.05;

20*log10(512/norm(f(:)-uexact1(:)))   
    
    
     [y nois_ma2] = acwmf(f*255);    
for i=1:3
    
u = tvinpaint(f,lambda,D,[],[],[],@tvregsimpleplot);

E = u(:) -f(:);
gamma = 1;

D = (abs(u-f)>thres);
D = double(D);

ker = [1, 1, 1, 1; 1, 1, 1, 1;1, 1, 1, 1;1,1,1,1];
Dt  = imfilter(D,ker,'same');
D   = D.*(Dt>=1);
D1  = D.*(Dt>=6);
D2  = D.*(Dt<6);

if i < 4
     [y nois_ma] = acwmf(u*255);
     nois_ma = 1-(1-nois_ma2).*(1-nois_ma);
     D2   = D2.*nois_ma;
end
D   = 1-(1-D2).*(1-D1);

thres = thres*0.8;

figure(5)
subplot(1,2,1)
imshow(D,[])
drawnow
subplot(1,2,2)
lambda= min(3*lambda,1000);
20*log10(512/norm(u(:)-uexact1(:))) 

end

[u, ~] = acwmf(u*255);
 20*log10(512/norm(u(:)/255-uexact1(:)))   
imshow(u/255,[])
toc
title('Inpainted');
