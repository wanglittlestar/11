function [v] = solvev(a,b)
% solve the followign op:
% min_{v} <v,a>
% s.t. 0<= v <= 1, 
% sum(v) >= b  <=>   -sum(v) <= -b
% b is an integer
% lb = ones(n,1)*0;
% ub = ones(n,1)*1;
% [v,fobj,exitflag] = linprog(a,-1*ones(1,n),-b,[],[],lb,ub);
[m1,n1]=size(a);
a = a(:);
n = m1 * n1;
[~,index]=sort(a,'ascend');
v = zeros(n,1);

v(index(1:b))=1;
 
v=reshape(v,m1,n1);



