function [x] = threadholding_l1_w(a,b)
% This program solves the following OP:
% min_{x} 0.5*x'*x + a'*x + <b,abs(x)>
% Here we assume b>=0

x = -sign(a) .* max(0,abs(a)-b);
