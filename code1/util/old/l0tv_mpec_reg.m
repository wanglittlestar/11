function [U] = l0tv_mpec_reg(B,O,Kmap,Ktmap,p,lambda,LargestEig,acc,B_Clean)

% min_{u}  lambda  || Dxu Dyu ||_{p,1} + ||o.*(Ku-b)||_0
% min_{u,v} lambda || Dxu Dyu ||_{p,1} + <e,e-v>,  s.t. <o.*|Au-b|,v> = 0
% min_{u,v,r,s,z} lambda || r s ||_{p,1} + <e,e-v>, s.t. <o.*|z|,v> = 0, Au-b = z, Dxu = r, Dyu = s
% L(u,v,r,s,z) =  || r s ||_{p,1}  + <e,e-v>
%                                  + <piz, Au-b-z> + 0.5 alpha |Au-b-z|_2^2
%                                  + <pir, Dxu-r> + 0.5 beta |Dxu-r|_2^2
%                                  + <pis, Dyu-s> + 0.5 beta |Dyu-s|_2^2,
%                                  + <piv, o.*|z|.*v> + 0.5 rho |o.*|z|.*v|_2^2

U = B;
dx = difX(U);
dy = difY(U);
R = dx;
S = dy;
V = ones(size(B));
Z = randn(size(B));

% Multipliers
piz = randn(size(B));
pir = randn(size(B));
pis = randn(size(B));
pio = randn(size(B));

gamma = 0.5*(1+sqrt(5)); %  golden ratio

alpha = 1;
rho   = 1;
beta  = 1;

for iter = 1:1e6,
    
    % Update Z
    cof_A = alpha * ones(size(B)) + rho * V.*V.*O;
    cof_B = - piz - alpha * (Kmap(U) - B);
    cof_C = pio.*V.*O;
    Z=threadholding_l1_w(cof_B./cof_A,cof_C./cof_A);
    
    % Update U
    g1 = Ktmap(piz) + alpha*Ktmap(Kmap(U)-B-Z);
    g3 = divX(-pir+beta*R)  - beta*divX(dx);
    g4 = divY(-pis+beta*S ) - beta*divY(dy);
    gradU =     g1 + g3 + g4;
    Lip = beta*4 + beta * 4 + alpha * LargestEig;
    U   = boxproj(U - gradU/Lip);    
    dx  = difX(U);    dy  = difY(U);% Update dx and dy whenever U has changed
    
    % Update RS
    R1 = - pir/beta - dx; S1 = - pis/beta - dy;
    [R,S] = threadholding_RS(R1,S1,lambda,beta,p);
    
    % Update V
    cof_A = rho * Z.*Z.*O;
    cof_B = pio.*abs(Z).*O - ones(size(B));
    V = - cof_B./cof_A;
    V(V>1)=1;
    V(V<0)=0;
    
    % Statistics
    r1 = fnorm(Kmap(U) - B-Z);
    r2 = fnorm(dx-R)+ fnorm(dy-S);
    r3 = fnorm(V.*abs(Z).*O);
    all = r1 + r2 + r3;
    if(iter>10&&all<acc),break;end
    PSNR   = snr(U,B_Clean);
    fprintf('iter:%d, nearness:(%.1e %.1e %.1e), pen:(%.1e %.1e %.1e), %f\n',iter,r1,r2,r3,alpha,beta,rho,PSNR);
    
    % Update Multipliers
    piz = piz + gamma*alpha*(Kmap(U) - B-Z);
    pir = pir + gamma*beta*(dx-R);
    pis = pis + gamma*beta*(dy-S);
    pio = pio + gamma*rho*(V.*abs(Z).*O);
    
    if(~mod(iter,10)),
        if(r1>r2 && r1>r3),
            alpha = alpha * 5;
        end
        if(r2>r1 && r2>r3),
            beta = beta *5;
        end
        if(r3>r1 && r3>r2),
            rho = rho * 5;
        end
    end
    
end

function [x] = boxproj(x)
x(x<0)=0;
x(x>1)=1;
