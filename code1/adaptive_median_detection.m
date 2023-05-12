function [X] = adaptive_median_detection(g1,g)

% This function finds the poisition of image corrupted by salt and pepper noise

% Input
%

% g1     --------------------   blur  image
% g       --------------------  blur and noise image
% lever --------------------   noise lever

% Ouput
% X ------  characteristic matrix 


%f=imread('C:\Program Files\MATLAB\R2013a\bin\work\图像复原\Penguins.jpg');
%f = imread('onion.png');
%  f = double(imread('House256.tif'))/255; 
%  H = fspecial('gaussian',9,10);
%  g1 = imfilter(f,H);
 image_gray = g1;
%image_gray=rgb2gray(f);%得到灰度图像
 ff = image_gray;
 ff(:) = 0;
alreadyProcessed = false(size(image_gray));%生成逻辑非的矩阵

% 迭代.
Smax=7;
for k = 3:2:Smax
   zmin = ordfilt2(image_gray, 1, ones(k, k), 'symmetric');
   zmax = ordfilt2(image_gray, k * k, ones(k, k), 'symmetric');
   zmed = medfilt2(image_gray, [k k], 'symmetric');
   
   processUsingLevelB = (zmed > zmin) & (zmax > zmed) & ...
       ~alreadyProcessed;
   zB = (image_gray > zmin) & (zmax > image_gray);
   outputZxy  = processUsingLevelB & zB;
   outputZmed = processUsingLevelB & ~zB;
   ff(outputZxy) = image_gray(outputZxy);
   ff(outputZmed) = zmed(outputZmed);
   
   alreadyProcessed = alreadyProcessed | processUsingLevelB;
   if all(alreadyProcessed(:))
      break;
   end
end

ff(~alreadyProcessed) = zmed(~alreadyProcessed);
f1 = g;


[m,n] = size(g); 
X = zeros(m,n);
for i = 1 : m
    for j = 1 : n
        if ff(i,j) ~= f1(i,j) & f1(i,j) == 0 || f1(i,j) == 1;
            X(i,j) = 0;
        else
            X(i,j) = 1;
        end
    end
end

            
