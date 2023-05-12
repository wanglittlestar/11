function test_quad_diag
n = 10;
a = rand(n,1);
b = randn(n,1);
k = 4;
x1=quadprog(diag(a),b,-ones(1,n),-k,[],[],zeros(n,1),ones(n,1))
sum(x1)

 
function [v] = quaddiag(a,b,k)
% solve the following OP:
% min_{v} 0.5 v'Av + b'v
% s.t. 0<=v<=1, v'e >= k
% where A is a diagonal matrix

% L(v,pil,piu,pi) =  0.5 v'Av + b'v + <pil,-v> + <piu,v-e> + pi (k-v'e)

