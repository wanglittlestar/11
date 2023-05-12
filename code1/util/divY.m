function R = divY(Py)
% return divergence matrix divP
% Px  : m x n
% divP: m x n  
R = Py-Py([1 1:end-1],:,:);
R(1,:,:)   = Py(1,:,:);        
R(end,:,:) = - Py(end-1,:,:);
