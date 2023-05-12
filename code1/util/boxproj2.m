function [x] = boxproj2(x)
x(x<0)=0; x(x>1)=1;

