function [U] = l0tv_padmm_color(B,O,Kmap,Ktmap,p,lambda,LargestEig,acc,B_Clean)

% min_{u}  lambda  sum_c=1^3 || Dxu Dyu ||_{p,1} + sum_c=1^3 ||o.*(Ku-b)||_0
% min_{u,v} lambda sum_c=1^3 || Dxu Dyu ||_{p,1} + sum_c=1^3 <e,e-v>,  s.t. <o.*|Au-b|,v> = 0
% min_{u,v,r,s,z} sum_c=1^3 lambda || r s ||_{p,1} + sum_c=1^3 <e,e-v>, s.t. <o.*|z|,v> = 0, Au-b = z, Dxu = r, Dyu = s
% L(u,v,r,s,z) =  sum_c=1^3 lambda || r s ||_{p,1}  + sum_c=1^3 <e,e-v>
%                                  + <piz, Au-b-z> + 0.5 alpha |Au-b-z|_2^2
%                                  + <pir, Dxu-r> + 0.5 beta |Dxu-r|_2^2
%                                  + <pis, Dyu-s> + 0.5 beta |Dyu-s|_2^2,
%                                  + <piv, o.*|z|.*v> + 0.5 rho |o.*|z|.*v|_2^2

[n1,n2,n3]=size(B);
E = ones(size(B));
B = double(B);
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
his = [];
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
    dx  = difX(U);    dy  = difY(U);% Update dx and dy whenever U has been changed
    
    % Update RS
    R1 = - pir/beta - dx; S1 = - pis/beta - dy;
    [R,S]=threadholding_RS(R1,S1,lambda,beta,p);
 
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
    
    D1 = Kmap(U) - B-Z;
    D21 = dx-R;
    D22 = dy-S;
    D3 = V.*abs(Z).*O;
    lobj = 0;
    for i=1:n3,
        lobj = lobj + lambda * sum(sum((R(:,:,i).^p + S(:,:,i).^p).^(1/p)));
        lobj = lobj + mdot(E,E-V(:,:,i));
        lobj = lobj + mdot(piz(:,:,i), D1(:,:,i)) + 0.5 * alpha * mdot(D1(:,:,i),D1(:,:,i));
        lobj = lobj + mdot(pir(:,:,i), D21(:,:,i)) + 0.5 * beta * mdot(D21(:,:,i),D21(:,:,i));
        lobj = lobj + mdot(pis(:,:,i), D22(:,:,i)) + 0.5 * beta * mdot(D22(:,:,i),D22(:,:,i));
        lobj = lobj + mdot(pio(:,:,i), D3(:,:,i)) + 0.5 * rho * mdot(D3(:,:,i),D3(:,:,i));
    end
    lobj
    his = [his;lobj];
    plot(his)
    pause(1)

    % Update Multipliers
%     piz = piz + gamma*alpha*(Kmap(U) - B-Z);
%     pir = pir + gamma*beta*(dx-R);
%     pis = pis + gamma*beta*(dy-S);
%     pio = pio + gamma*rho*(V.*abs(Z).*O);
%     
%     if(~mod(iter,10)),
%         if(r1>r2 && r1>r3),
%             alpha = alpha * 5;
%         end
%         if(r2>r1 && r2>r3),
%             beta = beta *5;
%         end
%         if(r3>r1 && r3>r2),
%             rho = rho * 5;
%         end
%     end
    
end

