function [obj] = maxnorm(a)
a = abs(a);
obj = max(a(:));