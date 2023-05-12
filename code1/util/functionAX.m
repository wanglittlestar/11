function [R] = functionAX(A,X,type)
if(strcmp(type,'denoising'))
    R = X;
elseif(strcmp(type,'deblurring'))
   R=imfilter(X,A,'circular');
%      R= conv2padded(X,A);
else
    error('unknown type');
end
