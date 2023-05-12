function [U] = l0tv_proj_cont(B,Amap,Atmap,p,kkk,LargestEig,acc,B_Clean)


% Primal Variables
U = B;

dx = difX(U);
dy = difY(U);
Z = 0*randn(size(B));
R = dx;
S = dy;

% Multipliers
piz = 0*randn(size(B));
pir = 0*randn(size(B));
pis = 0*randn(size(B));

gamma = 1.62; %  golden ratio
% gamma = 1.9;
p = 2;

alpha = 1000;
beta  = 1;

his = [];
for iter = 1:10000,
    
    % Update RS
    % min_{R,S} sum(sum((R.^p + S.^p).^(1/p))) + mdot(pir, -R) + 0.5*beta*fnorm(dx-R)^2 + mdot(pis, -S) + 0.5*beta*fnorm(dy-S)^2;
    % min_{R,S} sum(sum((R.^p + S.^p).^(1/p))) + 0.5*beta*fnorm(R-(dx+pir/beta))^2 + 0.5*beta*fnorm(S-(dy+pis/beta))^2;
    % min_{R,S} 1/beta*sum(sum((R.^p + S.^p).^(1/p))) + 0.5*fnorm(R-(dx+pir/beta))^2 + 0.5*fnorm(S-(dy+pis/beta))^2;
    R1 = - pir/beta -dx;    S1 = - pis/beta - dy;
    [R,S] = threadholding_RS(R1,S1,1,beta,p);
    
    
    % Update U
    dx  = difX(U);    dy  = difY(U);
    g1 = Atmap(piz) + alpha*Atmap(Amap(U)-B-Z);
    g2 = divX(-pir+beta*R)  - beta*divX(dx);
    g3 = divY(-pis+beta*S) - beta*divY(dy);
    gradU =     g1 + g2 + g3;
    Lip = beta*4 + beta*4 + alpha*LargestEig;
    U   = boxproj(U - gradU/Lip);
    dx  = difX(U);    dy  = difY(U);
    
    
    % Update Z:
    % min_{Z} mdot(piz, Amap(U)-B-Z) + 0.5 * alpha*fnorm(Amap(U)-B-Z)^2
    % min_{Z} mdot(piz, -Z/alpha) + 0.5 * fnorm(Z-(Amap(U)-B))^2
    % min_{Z} 0.5 * fnorm(Z-(Amap(U)-B+piz/alpha))^2
    Z1 = Amap(U) - B + piz/alpha;
    Z = proj_l0(Z1,kkk);
    
    
    Lag = sum(sum((R.^p + S.^p).^(1/p))) + mdot(piz, Amap(U)-B-Z) + 0.5 * alpha*fnorm(Amap(U)-B-Z)^2 ...
        + mdot(pir, dx-R) + 0.5*beta*fnorm(dx-R)^2 + mdot(pis, dy-S) + 0.5*beta*fnorm(dy-S)^2;
    his = [his;Lag];
    
    
    % Statistics
    r1 = fnorm(Amap(U) - B-Z);
    r2 = fnorm(dx-R)+ fnorm(dy-S);
    all = r1 + r2;
    
    if(iter>10&&all<0.1/255),break;end
    
    PSNR   =  snr(U, B_Clean);
    %     if(~mod(iter,100))
    fprintf('iter:%d, nearness:%.1e (%.1e %.1e), penalty:(%f %f), psnr: %f\n',iter,all,r1,r2,alpha,beta,PSNR);
    %     end
    
    % Update Multipliers
    piz = piz + gamma*alpha*(Amap(U) - B-Z);
    pir = pir + gamma*beta*(dx-R);
    pis = pis + gamma*beta*(dy-S);
    if(~mod(iter,30)),
        if(r1>r2),
            alpha = alpha * 2;
        else
            beta = beta * 2;
        end
    end
    
    
    
end
