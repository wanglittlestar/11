function [x] = threadholding_l0(y,lambda)
% solving the following OP:
% min_{x} 0.5 x'x + a'x + lambda * |x|_0
 

x = zeros(size(y));
I = find( ((y.^2)./(2*lambda))>1);
z = -y;
x(I) = z(I);
