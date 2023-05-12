function [obj] = fnorm(A)
[a,b,c]=size(A);

obj = 0;
for i=1:c,
obj = obj + norm(A(:,:,i),'fro');
end
