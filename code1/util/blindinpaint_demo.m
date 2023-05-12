% Blind inpaiting 
%
% This is a MATLAB script.  To run it, enter on the MATLAB console
%   >> blindinpaint_demo




lambda  = 2.1e1;
uexact  = phantom(512);
uexact  = uexact/max(uexact(:));
uexact1 = uexact;
uexact  = uexact+0.*randn(size(uexact));
uexact  = imnoise(uexact,'salt & pepper', 0.3);

% Construct inpainting region
[x,y]   = meshgrid(1:size(uexact,2),1:size(uexact,1));
[th,r]  = cart2pol(x-size(uexact,2)/2,y-size(uexact,1)/2);
D       = (sin(r/20.5+th) > 0.995);

f       = uexact;
f(D)    = rand(nnz(D),1);
%f       = f2;
f       = f3;
% Make a figure 
% clf;
% set(gcf,'Color',[1,1,1],'NumberTitle','off','Name','TV Inpainting');
% compareimages(f,'Input',f,'Inpainted');
% shg;

figure(2)
subplot(2,2,3)
imshow(f,[])
subplot(2,2,4)
imshow(uexact1,[])
[y nois_ma2] = acwmf(f*255);   
20*log10(512/norm(f(:)-uexact1(:)))

% Inpaint
D       = zeros(size(uexact,1),size(uexact,2));
thres   = 0.05;     
for i = 1:6
    u = tvinpaint(f,lambda,D,[],[],[],[]);
    E = u(:) -f(:);
    D = (abs(u-f)>thres);
    D = double(D);

% if i < 7
%      [y nois_ma] = acwmf(u*255);
%      nois_ma = 1-(1-nois_ma2).*(1-nois_ma);
%      D   = D.*nois_ma;
% end    
    
%    thres = thres*0.935;
    thres = thres*0.96;

%    figure(2)
    subplot(2,2,1)
    imshow(D,[])
    subplot(2,2,2)
    imshow(u,[])
    drawnow
%    lambda= min(1.27*lambda,1000);
    lambda= min(1.3*lambda,1000);
20*log10(512/norm(u(:)-uexact1(:)))
end
%title('Inpainted');
[u, ~] = acwmf(u*255);
20*log10(512/norm(u(:)/255-uexact1(:)))  
