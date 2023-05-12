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




function [U,total_iter,his,snrs] = l0tv_proj_cont2(B,Amap,Atmap,p,lambda,LargestEig,acc,B_Clean)

his = [];
snrs = [];

[m,n] = size(B);
V = ones(m,n);


% Primal Variables for the Subproblem
U = B;
[Dx] = difX(U);
[Dy] = difY(U);

% Dual Variables for the Subproblem
Z = ones(m,n);
R = difX(B);
S = difY(B);

% Multipliers
piz = zeros(m,n);
pir = zeros(m,n);
pis = zeros(m,n);

gamma = 1.6;
rho = 1;


total_iter = 0;

for iter=1:30,
    
    beta  = rho;
    alpha = rho*5;
    
    % Solve the folowing reweighted L1 subproblem via proximal ADMM
    % min_u   <e,e-v> + lambda sum(sum((Dx.^p + Dy.^p).^(1/p))) + rho <V,|Amap(U)-B|>;
    %         s.t. 0<=u<=1
    
    % Introducing some auxiliary equalities: R = Dx, S = Dy and Z = Amap(U)-B, we have the following
    % OP:
    % min_u   <e,e-v> + lambda sum(sum((R.^p + S.^p).^(1/p))) + rho <V,|Z|>;
    %         s.t. 0<=u<=1, R = Dx, S = Dy, Z = Amap(U)-B
    
    % The Lagrangian multiplier sequence is:
    % L(u,R,S,Z) = <e,e-v> + lambda sum(sum((R.^p + S.^p).^(1/p))) + rho
    % <V,|Z|> + <piz, Amap(U)-B> + 0.5 alpha*fnorm(Amap(U)-B-Z)^2
    % + <pir, Dx-R> + 0.5*beta*fnorm(Dx-R)^2;
    % + <pis, Dy-S> + 0.5*beta*fnorm(Dy-S)^2;
    
    for isub = 1:30,
        % Update Z
        C = Amap(U) - B;
        Z = threadholding_l1_w(-C-piz/alpha, (rho/alpha)*V );
        
        % Update RS
        R1 = - pir/beta -Dx;
        S1 = - pis/beta - Dy;
        if(p==2)
            normRS = sqrt(R1.^2 + S1.^2);normRS(normRS<1e-5)=1;
            Weight = (max(0,normRS-(lambda/beta))./ normRS);
            Weight(normRS==0)=1;
            R = -R1.*Weight;
            S = -S1.*Weight;
        else
            R = - sign(R1).*max(0,abs(R1)-(lambda/beta));
            S = - sign(S1).*max(0,abs(S1)-(lambda/beta));
        end
        
        % Update U
        last_U = U;
        g1 = Atmap(piz) + alpha*Atmap(Amap(U)-B-Z);
        g2 = divX(-pir+beta*R) - beta*divX(Dx);
        g3 = divY(-pis+beta*S) - beta*divY(Dy);
        gradU =     g1 + g2 + g3;
        total_iter = total_iter + 1 ;
        
        Lip = beta*4 + beta * 4 + alpha * LargestEig;
        U = min(max(U - gradU/Lip,0),1);
        
        his = [his;ComputeTrue(U,Amap,B,lambda,p)];
        sss = snr(U,B_Clean);
        snrs = [snrs;sss];
        
        [Dx] = difX(U);
        [Dy] = difY(U);
        
        stop = max(max(abs(last_U - U)));
        if(isub > 10 && stop < acc),break;end
        
        nearness1 = 0.5*(mnorm(Dx-R) + mnorm(Dy-S));
        nearness2 = mnorm(Amap(U) - B -Z);
        
        pir = pir + gamma*beta*(Dx-R);
        pis = pis + gamma*beta*(Dy-S);
        piz = piz + gamma*alpha*(Amap(U) - B -Z);
        
        
        
        
    end
    
    
    % Update V
    aaaa = max(max(V.*abs(Amap(U) - B)));
    PSNR   = snr(U,B_Clean);
    fprintf('iter:%d, rho: %d, constraint:%.1e, %f\n',iter,rho,aaaa,PSNR);
    if( aaaa<acc && nearness1<acc && nearness2<acc),
        break;
    end
    V = ones(m,n);
    V(abs(Amap(U) - B)*rho>1)=0;
    rho = rho * 2;
    
end


function [U] = l0tv0(B,Amap,Atmap,p,kkk,LargestEig,acc,B_Clean)


% Primal Variables
U = randn(size(B));

dx = difX(U);
dy = difY(U);
Z = randn(size(B));
R = dx;
S = dy;

% Multipliers
piz = 0.1*randn(size(B));
pir = 0.1*randn(size(B));
pis = 0.1*randn(size(B));

gamma = 1.62; %  golden ratio
gamma = 1.8;
p = 2;

alpha = 0.1;
beta  = 0.1;

his = [];
for iter = 1:10000,
    
    
    % Update Z:
    % min_{Z} mdot(piz, Amap(U)-B-Z) + 0.5 * alpha*fnorm(Amap(U)-B-Z)^2
    % min_{Z} mdot(piz, -Z/alpha) + 0.5 * fnorm(Z-(Amap(U)-B))^2
    % min_{Z} 0.5 * fnorm(Z-(Amap(U)-B+piz/alpha))^2
    Z1 = Amap(U) - B + piz/alpha;
    Z = proj_l1(Z1,kkk);
    
    
    % Update RS
    % min_{R,S} sum(sum((R.^p + S.^p).^(1/p))) + mdot(pir, -R) + 0.5*beta*fnorm(dx-R)^2 + mdot(pis, -S) + 0.5*beta*fnorm(dy-S)^2;
    % min_{R,S} sum(sum((R.^p + S.^p).^(1/p))) + 0.5*beta*fnorm(R-(dx+pir/beta))^2 + 0.5*beta*fnorm(S-(dy+pis/beta))^2;
    % min_{R,S} 1/beta*sum(sum((R.^p + S.^p).^(1/p))) + 0.5*fnorm(R-(dx+pir/beta))^2 + 0.5*fnorm(S-(dy+pis/beta))^2;
    
    R1 = - pir/beta -dx;
    S1 = - pis/beta - dy;
    if(p==2)
        normRS = sqrt(R1.^2 + S1.^2);normRS(normRS<1e-5)=1;
        Weight = (max(0,normRS-(1/beta))./ normRS);
        Weight(normRS==0)=1;
        R = -R1.*Weight;
        S = -S1.*Weight;
    else
        R = - sign(R1).*max(0,abs(R1)-(1/beta));
        S = - sign(S1).*max(0,abs(S1)-(1/beta));
    end
    
    
    
    
    % Update U
    dx  = difX(U);    dy  = difY(U);
    g1 = Atmap(piz) + alpha*Atmap(Amap(U)-B-Z);
    g2 = divX(-pir+beta*R)  - beta*divX(dx);
    g3 = divY(-pis+beta*S) - beta*divY(dy);
    gradU =     g1 + g2 + g3;
    Lip = beta*4 + beta*4 + alpha*LargestEig;
    U   = (U - gradU/Lip);
    dx  = difX(U);    dy  = difY(U);
    
    
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
    fprintf('iter:%d, nearness:%.1e (%.1e %.1e), %f %f %f\n',iter,all,r1,r2,PSNR,alpha,beta);
    %     end
    
    % Update Multipliers
    piz = piz + gamma*alpha*(Amap(U) - B-Z);
    pir = pir + gamma*beta*(dx-R);
    pis = pis + gamma*beta*(dy-S);
    if(~mod(iter,30)),
        if(r1>r2),
            alpha = alpha * 10;
        else
            beta = beta * 10;
        end
    end
    
    
    
end

function [U] = l0tv_mpec_cont(B,Kmap,Ktmap,p,kkk,LargestEig,acc,B_Clean)

% min_{u}   || Dxu Dyu ||_{p,1}
%      s.t. ||Ku-b||_0 <= k
% min_{u,v} || Dxu Dyu ||_{p,1}
%      s.t. <e,e-v> <= k,  <|Au-b|,v> = 0
% min_{u,v,r,s,z} || r s ||_{p,1}
%      s.t. <e,e-v> <= k,  <|z|,v> = 0, Au-b = z, Dxu = r, Dyu = s
% min_{u,v,r,s,z} || r s ||_{p,1}
%      s.t. <e,e-v> <= k,  <|z|,v> = 0, Au-b = z, Dxu = r, Dyu = s

% L(u,v,r,s,z) =  || r s ||_{p,1} + <piz, Au-b-z> + 0.5 alpha |Au-b-z|_2^2
%                                 + <pir, Dxu-r> + 0.5 beta |Dxu-r|_2^2
%                                 + <pis, Dyu-s> + 0.5 beta |Dyu-s|_2^2,
%                                 + <piv, |z|.*v> + rho <|z|,v>

% Primal Variables
[m,n]= size(B);
U = B;
V = ones(size(B));
dx = difX(U); dy = difY(U);
Z = 0.1*randn(size(B));
R = dx;
S = dy;

% Multipliers
piz = 0.0*randn(size(B));
pir = 0.0*randn(size(B));
pis = 0.0*randn(size(B));
piv = 0.0*rand(size(B));

gamma = 1.62; %  golden ratio
gamma = 1.8;
p = 2;

alpha = 1;
beta  = 1;
rho =  .1;

his = [];
for iter = 1:10000,
    
    for inner = 1:1,
        % Update Z:
        % o = Au-b;
        % min_{z} <piz, -z> + 0.5 alpha |z-o|_2^2 + <piv, |z|.*v> + rho <|z|,v>
        % min_{z} 0.5 alpha |z-(o+piz/alpha)|_2^2 + <piv.*v + rho v, |z|>
        % min_{z} 0.5 |z-(o+piz/alpha)|_2^2 + <piv.*v + rho v, |z|>/alpha
        o = Kmap(U) - B;
        o = o+piz/alpha;
        [Z] = threadholding_l1_w(-o,(piv.*V + rho*V)/alpha);
        
        
        % Update RS
        % min_{R,S} sum(sum((R.^p + S.^p).^(1/p))) + mdot(pir, -R) + 0.5*beta*fnorm(dx-R)^2 + mdot(pis, -S) + 0.5*beta*fnorm(dy-S)^2;
        % min_{R,S} sum(sum((R.^p + S.^p).^(1/p))) + 0.5*beta*fnorm(R-(dx+pir/beta))^2 + 0.5*beta*fnorm(S-(dy+pis/beta))^2;
        R1 = - pir/beta - dx;
        S1 = - pis/beta - dy;
        [R,S] = threadholding_RS(R1,S1,1,beta,p);
        
        % Update U
        % min_{u}     <piz, Au> + 0.5 alpha |Au-b-z|_2^2  + <pir, Dxu-r> + 0.5 beta |Dxu-r|_2^2  + <pis, Dyu-s> + 0.5 beta |Dyu-s|_2^2,
        
        dx  = difX(U);    dy  = difY(U);
        g1 = Ktmap(piz) + alpha*Ktmap(Kmap(U)-B-Z);
        g2 = divX(-pir+beta*R)  - beta*divX(dx);
        g3 = divY(-pis+beta*S) - beta*divY(dy);
        gradU =     g1 + g2 + g3;
        Lip = beta*4 + beta*4 + alpha*LargestEig;
        U   = boxproj(U - gradU/Lip);
        dx  = difX(U);    dy  = difY(U);
        
    end
    
    % Update V
    % min_{0<=v<=1}    <|z|.*piv + rho*|z|, v>   s.t.   <e,e-v> <= k       <=>
    V = solvev(abs(Z).*piv + rho*abs(Z),m*n-kkk);
    
    
    Lag = sum(sum((R.^p + S.^p).^(1/p))) + mdot(piz, Kmap(U)-B-Z) + 0.5 * alpha*fnorm(Kmap(U)-B-Z)^2 ...
        + mdot(pir, dx-R) + 0.5*beta*fnorm(dx-R)^2 + mdot(pis, dy-S) + 0.5*beta*fnorm(dy-S)^2 ...
        + mdot(piv,abs(Z).*V) + rho * mdot(abs(Z),V);
    his = [his;Lag];
    
    
    % Statistics
    r1 = fnorm(Kmap(U) - B-Z);
    r2 = fnorm(dx-R)+ fnorm(dy-S);
    r3 = mdot(abs(Z),V);
    all = r1 + r2 + r3;
    
    if(iter>10&&all<0.1/255),break;end
    PSNR   =  snr(U, B_Clean);
    %     if(~mod(iter,100))
    fprintf('iter:%d, nearness: (%f %f %f), penalty: (%f %f %f), psnr:%f\n',...
        iter,r1,r2,r3,alpha,beta,rho,PSNR);
    %     end
    
    
    %     plot(his)
    
    % Update Multipliers
    piz = piz + gamma*alpha*(Kmap(U) - B-Z);
    pir = pir + gamma*beta*(dx-R);
    pis = pis + gamma*beta*(dy-S);
    piv = piv + gamma*rho*(abs(Z).*V);
    
    if(~mod(iter,20)),
        imshow(U)
        pause(1)
        if(r1>r2 && r1>r3),
            alpha = alpha * 2;
        elseif(r2>r1 && r2>r3),
            beta = beta * 2;
        elseif(r3>r1 && r3>r2),
            rho = rho * 2;
        end
    end
    
    
    
end

function [U] = l0tv_mpec_cont2(B,Kmap,Ktmap,p,kkk,LargestEig,acc,B_Clean)

% min_{u}   || Dxu Dyu ||_{p,1}
%      s.t. ||Ku-b||_0 <= k
% min_{u,v} || Dxu Dyu ||_{p,1}
%      s.t. <e,e-v> <= k,  <|Au-b|,v> = 0
% min_{u,v,r,s,z} || r s ||_{p,1}
%      s.t. <e,e-v> <= k,  <|z|,v> = 0, Au-b = z, Dxu = r, Dyu = s
% min_{u,v,r,s,z} || r s ||_{p,1}
%      s.t. <e,e-v> <= k,  <|z|,v> = 0, Au-b = z, Dxu = r, Dyu = s

% L(u,v,r,s,z) =  || r s ||_{p,1} + <piz, Au-b-z> + 0.5 alpha |Au-b-z|_2^2
%                                 + <pir, Dxu-r>  + 0.5 beta  |Dxu-r|_2^2
%                                 + <pis, Dyu-s>  + 0.5 beta  |Dyu-s|_2^2,
%                                 + <piv, |z|.*v> + 0.5 rho   |z.*v|_2^2

% Primal Variables
[m,n]= size(B);
U = B;
V = ones(size(B));
dx = difX(U); dy = difY(U);
Z = 0.1*randn(size(B));
R = dx;
S = dy;

% Multipliers
piz = 0.0*randn(size(B));
pir = 0.0*randn(size(B));
pis = 0.0*randn(size(B));
piv = 0.0*rand(size(B));

gamma = 1.62; %  golden ratio
gamma = 1.8;
p = 2;

alpha = 1;
beta  = 1;
rho =  .1;

his = [];
for iter = 1:10000,
    
    % Update Z:
    % o = Au-b;
    % min_{z} <piz, -z> + 0.5 alpha |z-o|_2^2 + <piv, |z|.*v> + 0.5 rho |z.*v|_2^2
    % min_{z} 0.5 alpha |z-(o+piz/alpha)|_2^2 + <piv.*v,|z|> + 0.5 rho |z.*v|_2^2
    % y = o+piz/alpha, h=piv.*v
    % min_{z} 0.5 alpha |z-y|_2^2 + <h,|z|> + 0.5 rho |z.*v|_2^2
    o = Kmap(U) - B;
    y1 = o+piz/alpha;
    h1 = piv.*V;
    AAA = alpha  + ones(size(B)) + rho * V.*V;
    BBB = - alpha * y1;
    CCC = h1;
    BBB = BBB ./ AAA;
    CCC = CCC ./ AAA;
    [Z] = threadholding_l1(BBB,CCC);
    
    
    % Update RS
    % min_{R,S} sum(sum((R.^p + S.^p).^(1/p))) + mdot(pir, -R) + 0.5*beta*fnorm(dx-R)^2 + mdot(pis, -S) + 0.5*beta*fnorm(dy-S)^2;
    % min_{R,S} sum(sum((R.^p + S.^p).^(1/p))) + 0.5*beta*fnorm(R-(dx+pir/beta))^2 + 0.5*beta*fnorm(S-(dy+pis/beta))^2;
    R1 = - pir/beta - dx;
    S1 = - pis/beta - dy;
    [R,S] = threadholding_RS(R1,S1,1,beta,p);
    
    % Update U
    % min_{u}     <piz, Au> + 0.5 alpha |Au-b-z|_2^2  + <pir, Dxu-r> + 0.5 beta |Dxu-r|_2^2  + <pis, Dyu-s> + 0.5 beta |Dyu-s|_2^2,
    
    dx  = difX(U);    dy  = difY(U);
    g1 = Ktmap(piz) + alpha*Ktmap(Kmap(U)-B-Z);
    g2 = divX(-pir+beta*R)  - beta*divX(dx);
    g3 = divY(-pis+beta*S) - beta*divY(dy);
    gradU =     g1 + g2 + g3;
    Lip = beta*4 + beta*4 + alpha*LargestEig;
    U   = boxproj(U - gradU/Lip);
    dx  = difX(U);    dy  = difY(U);
    
    % Update V
    % min_{0<=v<=1}    <|z|.*piv + rho*|z|, v>   s.t.   <e,e-v> <= k       <=>
    V = solvev(abs(Z).*piv + rho*abs(Z),m*n-kkk);
    
    
    Lag = sum(sum((R.^p + S.^p).^(1/p))) + mdot(piz, Kmap(U)-B-Z) + 0.5 * alpha*fnorm(Kmap(U)-B-Z)^2 ...
        + mdot(pir, dx-R) + 0.5*beta*fnorm(dx-R)^2 + mdot(pis, dy-S) + 0.5*beta*fnorm(dy-S)^2 ...
        + mdot(piv,abs(Z).*V) + rho * mdot(abs(Z),V);
    his = [his;Lag];
    
    
    % Statistics
    r1 = fnorm(Kmap(U) - B-Z);
    r2 = fnorm(dx-R)+ fnorm(dy-S);
    r3 = mdot(abs(Z),V);
    all = r1 + r2 + r3;
    
    if(iter>10&&all<0.1/255),break;end
    PSNR   =  snr(U, B_Clean);
    %     if(~mod(iter,100))
    fprintf('iter:%d, nearness: (%f %f %f), penalty: (%f %f %f), psnr:%f\n',...
        iter,r1,r2,r3,alpha,beta,rho,PSNR);
    %     end
    
    
    %     plot(his)
    
    % Update Multipliers
    piz = piz + gamma*alpha*(Kmap(U) - B-Z);
    pir = pir + gamma*beta*(dx-R);
    pis = pis + gamma*beta*(dy-S);
    piv = piv + gamma*rho*(abs(Z).*V);
    
    if(~mod(iter,20)),
        imshow(U)
        pause(1)
        if(r1>r2 && r1>r3),
            alpha = alpha * 2;
        elseif(r2>r1 && r2>r3),
            beta = beta * 2;
        elseif(r3>r1 && r3>r2),
            rho = rho * 2;
        end
    end
    
    
    
end
