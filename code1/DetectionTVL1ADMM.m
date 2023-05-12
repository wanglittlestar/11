function [U,SNR,PSNR,SSIM,i,t2] = DetectionTVL1ADMM(I,F,H,O,maxit,espilon,mu,beta1,beta2,gamma)

%% % An Alternationg Direction Method of Multipliers (ADMM) is applied to solve the
% TVL1, which refers to the ADMM for solving TVL1 by J. yang, W. Yin, and Y. Zhang(2009)

%  consider L1 data fidelity term and TV regularization
% term,  The model is as follows:
%          min_{u} \|O.*(Ku-f)\|_1+\|Du\|_2      
%     O is the characteristic matrix
%=============================================================

%% Input
% H -------------------- convolution kernel representing K
% F -------------------- blurry and noise observation
% O -------------------- 0,1 matrix
% maxit ----------- represents the maximum number of  iteration steps 
% espilon ----------- Respectively represents the stop criteria 
% mu ---------------------- Regular parameters in the model
% beta1, beta2 -------------  Penalty parameter in ADMM
% gamma --------------------  Relaxation parameter in ADMM

%% Out 
% U ----------------------  Restore image data
% IPSNR, ISNR, ISSIM----------- Value of evaluation indicators after noise
% PSNR, SNR, SSIM------------ Value of evaluation indicators after recovery
% i ------------------- The number of steps required for iteration outside the algorithm
%% initialization
[m,n] = size(F);
C = getC(F,H); 
[D,Dt] = defDDt;
U = F;
Lam1 = zeros(m,n);
Lam2 = zeros(m,n);
Lam3 = zeros(m,n);
KUF = imfilter(U,H,'circular','conv') - F;
[D1U,D2U] = D(U);
Demon = beta1*C.eigsKtK+beta2*C.eigsDtD;
%% Algorithm ADMM
 tic;
t1 = clock;
for i = 1:maxit
    % X-subproblem
    X = prox_1(KUF+Lam1/beta1,mu/beta1,O);

    % Y-subproblem
    tv1 = D1U+Lam2/beta2;
    tv2 = D2U+Lam3/beta2;
    tvnorm = sqrt(tv1.^2+tv2.^2);
    tvnorm(tvnorm==0) = 1;
    tvnorm = max(tvnorm-1/beta2,0)./tvnorm;
    Y1 = tv1.*tvnorm;
    Y2 = tv2.*tvnorm;

    % U-subproblem
    Up = U;
    TKXL = beta1*X-Lam1;
    TKXL = imfilter(TKXL,H,'circular','corr');
    TKF = beta1*F;
    TKF = imfilter(TKF,H,'circular','corr');
    TDYL = Dt(beta2*Y1-Lam2,beta2*Y2-Lam3);
    U = TKXL+TKF+TDYL;
    U = fft2(U)./Demon;
    U = real(ifft2(U));

    % update Lam
    KUF = imfilter(U,H,'circular','conv')-F;
    [D1U,D2U] = D(U);
    Lam1 = Lam1+gamma*beta1*(KUF-X);
    Lam2 = Lam2+gamma*beta2*(D1U-Y1);
    Lam3 = Lam3+gamma*beta2*(D2U-Y2);

    SNR = 10*(log10(size(I,1)*size(I,2))-log10(norm(I(:)-U(:))^2));
    [PSNR(i),ppSNR]=psnr(I,U);
    t2(i)=etime(clock,t1); 
%    SNR(i) = 20*log10(norm(I(:))/norm(I(:)-U(:)));
    SSIM = ssim(I,U);
     fprintf('i:%d, Psnr:%0.4f,ssim:%0.4f\n',i,PSNR(i),SSIM);
%    f = fval(D1U,D2U,KUF);
     r1 = fnorm(KUF-X);
     r2 = fnorm(D1U-Y1)+ fnorm(D2U-Y2);
    if r1+r2<=espilon
        break
    end
end

%%
function C = getC(F,H)
    [m,n] = size(F);
    C.eigsDtD = abs(psf2otf([1,-1],[m,n])).^2 + abs(psf2otf([1;-1],[m,n])).^2;
    C.eigsK = psf2otf(H,[m,n]);
    C.eigsKt = conj(C.eigsK);
    C.eigsKtK = abs(C.eigsK).^2;
    C.proI = abs(psf2otf([1],[m,n]));
    C.KtF = real(ifft2(C.eigsKt .* fft2(F)));
return;

function [D,Dt] = defDDt
    D = @(U) ForwardD(U);
    Dt = @(X,Y) Dive(X,Y);

 return;

 function [Dux,Duy] = ForwardD(U)
        Dux = [diff(U,1,2), U(:,1) - U(:,end)];
        Duy = [diff(U,1,1); U(1,:) - U(end,:)];
 return;


 function DtXY = Dive(X,Y)
        DtXY = [X(:,end) - X(:, 1), -diff(X,1,2)];
        DtXY = DtXY + [Y(end,:) - Y(1, :); -diff(Y,1,1)];
 return;
 
%  function [f] = fval(D1U,D2U,KUF)
%        tv = sum(sum(sqrt(D1U.^2 + D2U.^2)));
%        fid = sum(sum(abs(KUF)));
%        f = tv + fid;
%  return;
%  

function X = prox_1(KUF,lammad,O)
[m,n] = size(KUF);
X =zeros(m,n);
for i = 1:m
    for j = 1:n
        if O(i,j) == 1
             X(i,j) = max(abs(KUF(i,j)) - lammad,0).*sign(KUF(i,j)); 
        else
            X(i,j) = KUF(i,j);
        end
    end
end
return

