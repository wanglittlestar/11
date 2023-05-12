function [U,peaksnr,SNR,iter,t2] = l0tv_padmm_color_for_plot_fobj11(B,O,Kmap,Ktmap,p,lambda,LargestEig,acc,B_Clean,pen_ratio,opts,miu)
sizeB = size(B);
B = double(B);
[peaksnr, snr] = psnr(B_Clean,B);
fprintf('\n The Peak-SNR value is %0.4f', peaksnr);
[ssimval, ssimmap] = ssim(B_Clean,B);
fprintf('The SSIM value is %0.4f.\n',ssimval);
U = B;
dx = difX(U);
dy = difY(U);
Kub = Kmap(U) - B;
V = ones(sizeB);
z1=ones(sizeB);

% Multipliers
piz = 1e-3*randn(sizeB);
pir = 1e-3*randn(sizeB);       
pis = 1e-3*randn(sizeB);
pio = 1e-3*randn(sizeB);
pim= 1e-3*randn(sizeB);
gamma = 0.5*(1+sqrt(5)); %  golden ratio
opts.alpha ;
opts.beta  ;
opts.rho1  ;
opts.rho2;
opts.iter;
ratio = pen_ratio;
alpha=opts.alpha;
beta=opts.beta;
rho1=opts.rho1;
rho2=opts.rho2;
iter=opts.iter;
 tic;
t1 = clock;
for iter= 1:iter
    % Update Z(y subproblem)
    VO = V.*O;
    cof_A = rho1 * V.*VO + alpha;
    cof_B = - piz - alpha * Kub;
    cof_C = pio.*VO;
    Z=threadholding_l1_w(cof_B./cof_A,cof_C./cof_A);
    
    % Update RS(x,subproblem)
    [R,S]=threadholding_RS(- pir/beta - dx,- pis/beta - dy,lambda,beta,p);
    
    % Update U
    g1 = Ktmap(piz) + alpha*Ktmap(Kub-Z);
    g3 = divX(-pir+beta*R)  - beta*divX(dx);
    g4 = divY(-pis+beta*S ) - beta*divY(dy);
    g5=pim+rho2*(U-z1);
    gradU =     g1 + g3 + g4+g5;    
    Lip = beta*4 + beta * 4 + alpha * LargestEig+rho2;
    U   = boxproj(U - gradU/Lip);
    
    dx  = difX(U);    dy  = difY(U);    Kub = Kmap(U) - B;
    % Update {dx,dy,Kub} whenever U has changed
    
    % Update V
    ZO = Z.*O;
    absZO = abs(ZO);
    cof_A = rho1 * ZO.*ZO;%%%%%
    cof_B = pio.*absZO - 1;
    V = boxproj(-cof_B./cof_A);
    %update z1
    A=U+pim./rho2;
    lambda1=miu/rho2;
    z1=Prox_lambda_nuclear_norm(A,lambda1);

    % Update Multipliers
    Kubz = Kub-Z;
    dxR = dx-R;
    dyS = dy-S;
    VabsZO = V.*absZO;    
    
    piz = piz + gamma*alpha*Kubz;
    pir = pir + gamma*beta*dxR;
    pis = pis + gamma*beta*dyS;
    pio = pio + gamma*rho1*VabsZO;
    pim = pim + gamma*rho2*(U-z1); 
    
    [peaksnr(iter),snr] = psnr(B_Clean,U);
      t2(iter)=etime(clock,t1); 
    [ssimval,ssimap] = ssim(B_Clean,U); 
    SNR=psnr(B_Clean,U);
    r1 = fnorm(Kubz);
    r2 = fnorm(dxR)+ fnorm(dyS);
    r3 = fnorm(VabsZO);
    r4=fnorm(U-z1);
    all = r1 + r2 + r3+r4;
    if(iter>30&&all<acc)
        break;
    end
    fprintf('iter:%d, Psnr:%0.4f,ssim:%0.4f\n',iter,peaksnr(iter),ssimval);
%     [peaksnr,snr] = psnr(B_Clean,U);
%     [ssimval,ssimap] = ssim(B_Clean,U); 
%     if(~mod(iter,30)),
%          fprintf('iter:%d, Psnr:%0.4f,ssim:%0.4f\n',iter,peaksnr,ssimval);
%     end
    
    if(~mod(iter,30)),
        if(r1>r2 && r1>r3&&r1>r4),
            alpha = alpha * ratio;
        end
        if(r2>r1 && r2>r3&&r2>r4),
            beta = beta * ratio;
        end
        if(r3>r1 && r3>r2&&r3>r4),
            rho1 = rho1 * ratio;
        end
        if(r4>r1 && r4>r2&&r4>r3),
            rho2 = rho2 * ratio;
        end
    end

end

function [R,S] = threadholding_RS(A,B,lambda,beta,p)
% Solve the following OP:
% min_{R,S} lambda sum(sum((R.^p + S.^p).^(1/p))) + 0.5 beta ||R + A||_F^2 + 0.5 beta ||S + B||_F^2

if(p==1)
     R = - sign(A).*max(0,abs(A)-(lambda/beta));
     S = - sign(B).*max(0,abs(B)-(lambda/beta));
elseif(p==2)
     normRS = sqrt(A.^2 + B.^2);
     Weight = max(0,1 - (lambda/beta)./ normRS);
     R = -A.*Weight;
     S = -B.*Weight; 
end
function [x] = threadholding_l1_w(a,b)
% This program solves the following OP:
% min_{x} 0.5*x'*x + a'*x + <b,abs(x)>
% Here we assume b>=0

x = -sign(a) .* max(0,abs(a)-b);
function [x] = boxproj(x)
x=max(x,0); x=min(x,1);
function [ Prox_lambda_nuclear_norm ] = Prox_lambda_nuclear_norm(A,lambda)
% iter=10001;
% i=1;
% A = rand(2,3);
% while i < iter
% % A = [1,2,3;
% %     -1,-2,-3];
%lambda = 1;


[s,u,v] = svd(A);
S_lambda = sign(u).*max(abs(u) - lambda,0);
Prox_lambda_nuclear_norm = s*S_lambda*v';
