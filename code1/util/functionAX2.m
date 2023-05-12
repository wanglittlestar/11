function [R] = functionAX(A,X,type)
if(strcmp(type,'denoising'))
    R = X;
elseif(strcmp(type,'deblurring'))
    R=imfilter(X,A,'symmetric');
else
    error('unknown type');
end
