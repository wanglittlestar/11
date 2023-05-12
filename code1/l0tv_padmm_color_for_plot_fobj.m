function [U,peaksnr,SNR,iter,t] = l0tv_padmm_color_for_plot_fobj(B,O,Kmap,Ktmap,p,lambda,LargestEig,acc,B_Clean,pen_ratio,opts)
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

% Multipliers
piz = 1e-3*randn(sizeB);
pir = 1e-3*randn(sizeB);       
pis = 1e-3*randn(sizeB);
pio = 1e-3*randn(sizeB);

gamma = 0.5*(1+sqrt(5)); %  golden ratio
opts.alpha ;
opts.beta  ;
opts.rho  ;
ratio = pen_ratio;
alpha=opts.alpha;
beta=opts.beta;
rho=opts.rho;
iter=opts.iter;
 tic;
t1 = clock;
for iter = 1:iter,
    
    % Update Z
    VO = V.*O;
    cof_A = rho * V.*VO + alpha;
    cof_B = - piz - alpha * Kub;
    cof_C = pio.*VO;
    Z=threadholding_l1_w(cof_B./cof_A,cof_C./cof_A);
    
    % Update RS
    [R,S]=threadholding_RS(- pir/beta - dx,- pis/beta - dy,lambda,beta,p);
    
    % Update U
    g1 = Ktmap(piz) + alpha*Ktmap(Kub-Z);
    g3 = divX(-pir+beta*R)  - beta*divX(dx);
    g4 = divY(-pis+beta*S ) - beta*divY(dy);
    gradU =     g1 + g3 + g4;
    Lip = beta*4 + beta * 4 + alpha * LargestEig;
    U   = boxproj(U - gradU/Lip);
    
    dx  = difX(U);    dy  = difY(U);    Kub = Kmap(U) - B;
    % Update {dx,dy,Kub} whenever U has changed
    
    % Update V
    ZO = Z.*O;
    absZO = abs(ZO);
    cof_A = rho * Z.*ZO+1;
    cof_B = pio.*absZO - 1-V;
    V = boxproj(-cof_B./cof_A);
    
    % Update Multipliers
    Kubz = Kub-Z;
    dxR = dx-R;
    dyS = dy-S;
    VabsZO = V.*absZO;    
    
    piz = piz + gamma*alpha*Kubz;
    pir = pir + gamma*beta*dxR;
    pis = pis + gamma*beta*dyS;
    pio = pio + gamma*rho*VabsZO;
   
    [peaksnr(iter), snr] = psnr(B_Clean,U);
     t(iter)=etime(clock,t1); 
    [ssimval, ssimmap] = ssim(B_Clean,U);
     SNR=psnr(B_Clean,U);
    % Statistics
    r1 = fnorm(Kubz);
    r2 = fnorm(dxR)+ fnorm(dyS);
    r3 = fnorm(VabsZO);
    all = r1 + r2 + r3;
    if(iter>30&&all<acc),break;end
     fprintf('iter:%d, Psnr:%0.4f,ssim:%0.4f\n',iter,peaksnr(iter),ssimval);
    
    if(~mod(iter,30)),
        if(r1>r2 && r1>r3),
            alpha = alpha * ratio;
        end
        if(r2>r1 && r2>r3),
            beta = beta * ratio;
        end
        if(r3>r1 && r3>r2),
            rho = rho * ratio;
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