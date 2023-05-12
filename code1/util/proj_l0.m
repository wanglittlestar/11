function [X] =  projl0(Y,k)
% solving the following OP:
% min_{x} 0.5 |x-y|_2^2
% s.t. |x|_0 <= k 

[m,n]=size(Y);
[x] =  projl0_vector(Y(:),k);
X = reshape(x,m,n);




function [x] =  projl0_vector(y,k)
% solving the following OP:
% min_{x} 0.5 |x-y|_2^2
% s.t. |x|_0 <= k 
y1 = abs(y);
[value,index]=sort(y1,'descend');
index = index(1:k);
x=zeros(size(y));
x(index) = y(index);