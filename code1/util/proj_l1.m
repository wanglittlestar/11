function  X = proj_l1(Z,alfa)
[m,n]=size(Z);
x = proj_l1_vector(Z(:),alfa);
X = reshape(x,m,n);


function w = proj_l1_vector(v, b)
%    min   ||w - v||_2
%    s.t.  ||w||_1 <= b.
try
u = sort(abs(v),'descend');
sv = cumsum(u);
rho = find(u > (sv - b) ./ (1:length(u))', 1, 'last');
theta = max(0, (sv(rho) - b) / rho);
% force w to a small number 1e-8 instead of strict zero
w = sign(v) .* max(abs(v) - theta,1e-8);
catch exception
    w = randn(size(v));
 
end