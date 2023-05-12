function [r] = nnz2(x)
x = x(:);
x = abs(x);
x=sort(x,'descend');
sum1 = 0;
epsi = 0.00001;
% epsi = 1/255;
% epsi = 0.001;
norm1 = (1-epsi)*sum(x);
n = length(x);

for i=1:n,
sum1 = sum1+x(i);
if(sum1 >=norm1),break;end
end
r = i;


    