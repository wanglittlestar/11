% TV inpainting demo using tvrestore, Pascal Getreuer 2010
%
% This is a MATLAB script.  To run it, enter on the MATLAB console
%   >> tvinpaint_demo

path(path,'C:\Users\yanm\Dropbox\BlindInpaiting\DECOLOR_MATLAB_WINDOWS\DECOLOR_MATLAB_WINDOWS\gco-v3.0\matlab')
path(path,'C:\Users\yanm\Dropbox\BlindInpaiting\DECOLOR_MATLAB_WINDOWS\DECOLOR_MATLAB_WINDOWS\internal')

%lambda = 1.25e1;
lambda = 1.25e1;
%lambda = 1.25;
uexact = double(sum(imread('boat.512.tiff'),3))/255;

uExact = double(sum(imread('lena512.bmp'),3));

uexact = uExact/255;

%uexact = phantom(512);
% D1 = (uexact == 1);
% uexact(D1) = 0;
%uexact = uexact/max(uexact(:));
%uexact = phantom(512);
%K = fspecial('gaussian',7,1);

uexact1 = uexact;
%uexact = conv2(uexact,K,'same');
fuexact = uexact+5/255*randn(size(uexact));


%uexact = imnoise(uexact,'salt & pepper', 0.3);
uexact   = addnoise(uexact*255, 40/100, 'rd')/255;

% Construct inpainting region
[x,y] = meshgrid(1:size(uexact,2),1:size(uexact,1));
[th,r] = cart2pol(x-size(uexact,2)/2,y-size(uexact,1)/2);
%D = (sin(r/1.5+th) > 0.950);

%D = (sin(r/20.5+th) > 0.993);
%D = (sin(r/20.5+th) > 0.995);
%D = (sin(r/20.5+th) > 0.997);
%D = (sin(r/20.5+th) > 0.999);

D = (sin(r/4.5+th) > 0.9);


f = uexact;

tic
%f(D) = 1*rand(nnz(D),1);

 PSNR  = 20*log10(512/norm(f(:)-uexact1(:)))
clear f3;

f = f3;
% Make a figure 
clf;
set(gcf,'Color',[1,1,1],'NumberTitle','off','Name','TV Inpainting');
compareimages(f,'Input',f,'Inpainted');
shg;

% Inpaint
D = zeros(size(uexact,1),size(uexact,2));
D2 =D;
    thres = 0.09;
 20*log10(512/norm(f(:)-uexact1(:)))   
    
%     hMRF = GCO_Create(numel(uexact),2);
%     GCO_SetSmoothCost( hMRF, [0 1;1 0] );
%         AdjMatrix = getAdj(size(uexact));
%     amplify = 10 * 3;
%     GCO_SetNeighbors( hMRF, amplify * AdjMatrix );
    
     [y nois_ma2] = acwmf(f*255);    
for i=1:2
    
u = tvinpaint(f,lambda,D,[],[],[],@tvregsimpleplot);
%u = tvinpaint(f,lambda,D,'l1',[],[],@tvregsimpleplot);

E = u(:) -f(:);
gamma = 1;

%         % call GCO to run graph cuts
%         energy_cut = 0;
%             GCO_SetDataCost( hMRF, (amplify/gamma)*[ 0.5*E.^2, ~D(:)*lambda + D(:)*0.5*max(E(:)).^2]' );
%             GCO_Expansion(hMRF);
%             D = ( GCO_GetLabeling(hMRF) ~= 1 )';


D = (abs(u-f)>thres);
D = double(D);

ker = [1, 1, 1, 1; 1, 1, 1, 1;1, 1, 1, 1;1,1,1,1];
Dt  = imfilter(D,ker,'same');
D   = D.*(Dt>=1);

% if i < 4
%      [y nois_ma] = acwmf(u*255);
%      nois_ma = 1-(1-nois_ma2).*(1-nois_ma);
%      D   = D.*nois_ma;
% end
%      

D2 = D;
     % if i < 5
%     D = D.*(Dt>=5);
% elseif i < 10
%     D = D.*(Dt>=5);
% elseif i <15
%     D = D.*(Dt>=4);
% elseif i < 20
%     D = D.*(Dt>=4);
% else
%     D = D.*(Dt>=3);
% end
        
% Dt  = imfilter(D,ker,'same');
% 
% D = D.*(Dt>=3);

thres = max(thres*0.4,0.007);


% if i == 10
%     thres = 0.03;
% elseif i == 15
%     thres = 0.03;
% elseif i == 20
%     thres = 0.03;
% end
figure(5)
subplot(1,2,1)
imshow(D,[])
drawnow
subplot(1,2,2)
lambda= min(2.8*lambda,1000);
 20*log10(512/norm(u(:)-uexact1(:))) 

% if i == 5
%     lambda = 4e1;
% elseif i == 10
%     lambda = 5e1;
% elseif i == 15
%     lambda = 6e1;
% elseif i == 20
%     lambda = 6e1;
% end
end

[u, ~] = acwmf(u*255);
 20*log10(512/norm(u(:)/255-uexact1(:)))   
imshow(u/255,[])
toc
title('Inpainted');
