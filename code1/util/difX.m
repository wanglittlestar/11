function [R] = difX(X)
R = X(:,[2:end end],:)-X;
% 在x方向进行差分