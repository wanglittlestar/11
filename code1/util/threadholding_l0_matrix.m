function [X] = threadholding_l0_matrix(Y,lambda)
% solving the following OP:
% min_{x} 0.5 x'x + a'x + lambda * |x|_0
 
X= threadholding_l0(Y(:),lambda(:));
X = reshape(X,size(Y));
