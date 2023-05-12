
function [out,peaksnr,SNR,k,t2]= L0_OGS_ADMM1(f,Img,K,opts)
% f=double(f);
[peaksnr, snr] = psnr(Img,f);
fprintf('\n The Peak-SNR value is %0.4f', peaksnr);
[ssimval, ssimmap] = ssim(Img,f);
fprintf('The SSIM value is %0.4f.\n',ssimval);

lam         = opts.lam; 
o           = opts.O; % mask
Nit         = opts.Nit;
tol         = opts.tol; 
grpSz       = opts.grpSz; %Group size
Nit_inner   = opts.Nit_inner;
beta_1 =10;
beta_2 = 1000;
beta_3 = 50000;
[row, col]  = size(f);
u           = f;
v = ones(size(u));
relError        = zeros(Nit,1);
psnrGain        = relError;     % PSNR improvement every iteration
ssimGain        = relError;

epsi_1         = 1e-3*randn(size(f));
epsi_2        = epsi_1;
delta         = epsi_1;
pio          = epsi_1;


eigK        = psf2otf(K,[row col]); %In the fourier domain
eigKtK      = abs(eigK).^2;
eigDtD      = abs(fft2([1 -1], row, col)).^2 + abs(fft2([1 -1]', row, col)).^2;

[D,Dt]      = defDDt(); %Declare forward finite difference operators
[Dux, Duy] = D(u);


lhs         = beta_2*eigKtK + beta_1*eigDtD; % From normal eqns.
Ku__f       = imfilter (u,K,'circular') -f; % Ku-f
 tic;
t1 = clock;
for k = 1:Nit
    
     u_old   = u;
     %*** solve y - subproblem ***
     w = v.*o;
      
     A = beta_3*v.*w + beta_2;
     B = -delta - beta_2*Ku__f;
     C = pio.*w;
     
     y = threshold_1(B./A,C./A);
     
     %*** solve x - subprob1lem (OGSTV problem) ***
     
     
     x1 = gstvdm(Dux + epsi_1/beta_1 , grpSz , lam/beta_1, Nit_inner);
     x2 = gstvdm(Duy + epsi_2/beta_1 , grpSz , lam/beta_1, Nit_inner);
      
    
     %*** solve u - subproblem ***
     ftemp   =  beta_2*y + beta_2*f - delta;
     rhs     = imfilter(ftemp,K,'circular') + Dt(beta_1*x1 - epsi_1 ,beta_1*x2 - epsi_2);
     u       = fft2(rhs)./lhs;
     u       = real(ifft2(u));
      
      
     [Dux, Duy]  = D(u);
     Ku__f       = imfilter (u,K,'circular') -f;
      
     %*** solve v - subproblem ***
     s = beta_3*y.*y.*o  ;
     
     c = pio.*abs(y.*o) - 1 ;
     v = max(-c./s,0);
     v = min(v,1);
         
     %*** Update the Lagrange multipliers ***
     epsi_1 = epsi_1 + beta_1*(Dux - x1);
     epsi_2 = epsi_2 + beta_1*(Duy - x2);
      
     delta = delta + beta_2*(Ku__f - y);
     pio   = pio   + beta_3*(v.*abs(y.*o));
      
  
  
     [peaksnr(k), snr] = psnr(Img,u);
      t2(k)=etime(clock,t1); 
     [ssimval, ssimmap] = ssim(Img,u); 
     SNR = 10*(log10(size(Img,1)*size(Img,2))-log10(norm(Img(:)-f(:))^2));
     fprintf('Nit:%d, Psnr:%0.4f,ssim:%0.4f\n',k,peaksnr(k),ssimval);
      r1 = fnorm(Ku__f - y);
     r2 = fnorm(Dux - x1)+ fnorm(Duy - x2);
     r3 = fnorm(v.*abs(y.*o));
     all = r1 + r2 + r3;
     if(Nit>30&&all<tol),break;end
end
    out= u;

end

function [D,Dt] = defDDt()
D  = @(U) ForwardDiff(U);
Dt = @(X,Y) Dive(X,Y);
end

function [Dux,Duy] = ForwardDiff(U)
 Dux = [diff(U,1,2), U(:,1,:) - U(:,end,:)];
 Duy = [diff(U,1,1); U(1,:,:) - U(end,:,:)];
end

function DtXY = Dive(X,Y)
  % Transpose of the forward finite difference operator
  % is the divergence fo the forward finite difference operator
  DtXY = [X(:,end) - X(:, 1), -diff(X,1,2)];
  DtXY = DtXY + [Y(end,:) - Y(1, :); -diff(Y,1,1)];   
end


function x = threshold_1(a,b)
x = -sign(a).*max(0,abs(a) - b);
end

