function [fobj] = ComputeTrue(U,Amap,B,lambda,p)
[R] = difX(U);
[S] = difY(U);
fobj = nnz2(Amap(U)-B) + lambda * sum(sum((R.^p + S.^p).^(1/p))) ;


% function [r] = nnz3(x)
% x = abs(x);
% the = 1/255;
% x(x<=the) = 0;
% x(x>the) = 0;
% r = sum(x);

