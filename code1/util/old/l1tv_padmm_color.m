function [U] = l1tv_padmm_color(B,O,Kmap,Ktmap,p,lambda,LargestEig,acc,B_Clean)
% min_{u}  lambda  sum_c=1^3 || Dxu Dyu ||_{p,1} + sum_c=1^3 ||o.*(Ku-b)||_1
% min_{u,r,s,z} sum_c=1^3 lambda || r s ||_{p,1} + sum_c=1^3 <o,|z|>, s.t. Au-b = z, Dxu = r, Dyu = s
% L(u,r,s,z) =  sum_c=1^3 lambda || r s ||_{p,1} + sum_c=1^3 <o,|z|>
%                                  + <piz, Au-b-z> + 0.5 alpha |Au-b-z|_2^2
%                                  + <pir, Dxu-r> + 0.5 beta |Dxu-r|_2^2
%                                  + <pis, Dyu-s> + 0.5 beta |Dyu-s|_2^2,

sizeB = size(B);
B = double(B);
U = B;
Z = B;
dx = difX(U);
dy = difY(U);
Kub = Kmap(U) - B;
R = zeros(size(B));
S = zeros(size(B));

const = 1;
% Multipliers
piz = const*randn(sizeB);
pir = const*randn(sizeB);
pis = const*randn(sizeB);

gamma = 1.6; %  golden ratio
alpha = 50;
beta  = 10;

snrs =[];
his=[];

for iter = 1:50,
    
    % Update Z
    % <o,|z|> + <piz, -z> + 0.5 alpha |Au-b-z|_2^2
    % <o,|z|> + <piz, -z> + 0.5 alpha |z-(Au-b)|_2^2
    % <o,|z|> + 0.5 alpha |z-(Au-b)-piz/alpha|_2^2
    % <o,|z|> + 0.5 alpha |z-(Au-b+piz/alpha)|_2^2
    Z=threadholding_l1_w(-(Kub+piz/alpha),O);
    [fobj] = computeObj(U,O,Kmap,p,Z,B,R,S,piz,pir,pis,lambda,alpha,beta);
    his = [his;fobj];
     plot(his)
   checkdec(his);
    
    % Update RS
    % lambda|| R S ||_{p,1} + <pir, -r> + 0.5 beta |r-Dxu|_2^2 + <pis, -s> + 0.5 beta |s-Dyu|_2^2,
    % lambda|| R S ||_{p,1} + 0.5 beta |r-Dxu-pir/beta|_2^2 + 0.5 beta |s-Dyu-pis/beta|_2^2,
    [R,S]=threadholding_RS(- pir/beta - dx,- pis/beta - dy,lambda,beta,p);
    [fobj] = computeObj(U,O,Kmap,p,Z,B,R,S,piz,pir,pis,lambda,alpha,beta);
    his = [his;fobj];
    plot(his)
    checkdec(his);
    
    % Update U
    %    <piz, Au> + 0.5 alpha |Au-b-z|_2^2 
    %  + <pir, Dxu-r> + 0.5 beta |Dxu-r|_2^2
    %  + <pis, Dyu-s> + 0.5 beta |Dyu-s|_2^2,

    g1 = Ktmap(piz) + alpha*Ktmap(Kub-Z);
    g3 = divX(-pir+beta*R)  - beta*divX(dx);
    g4 = divY(-pis+beta*S ) - beta*divY(dy);
    gradU =     g1 + g3 + g4;
    Lip = 100000;% beta*4 + beta * 4 + alpha * LargestEig;
    U   = boxproj(U - gradU/Lip);
    dx  = difX(U);    dy  = difY(U);    Kub = Kmap(U) - B;
    % Update {dx,dy,Kub} whenever U has changed
    
    [fobj] = computeObj(U,O,Kmap,p,Z,B,R,S,piz,pir,pis,lambda,alpha,beta);
    his = [his;fobj];
      plot(his)
  checkdec(his);


    % Update Multipliers
%      Kubz = Kub-Z;
%      dxR = dx-R;
%      dyS = dy-S;
% 
%      piz = piz + gamma*alpha*Kubz;
%      pir = pir + gamma*beta*dxR;
%      pis = pis + gamma*beta*dyS;

    % Statistics
%     r1 = fnorm(Kubz);
%     r2 = fnorm(dxR)+ fnorm(dyS);
%     all = r1 + r2;
%     if(iter>30&&all<acc),break;end
% fprintf('%f ',all);
%      PSNR   = snr_l1(U,B_Clean);
%      snrs = [snrs;PSNR];
%      his = [his;computeTrueObj(U,Kmap,B,O,p,lambda)];
%     
%     if(~mod(iter,30)),
%         fprintf('iter:%d, nearness:(%.1e %.1e %.1e), pen:(%.1e %.1e), %f\n',iter,r1,r2,all,alpha,beta,PSNR);
%     end
%     [fobj] = computeTrueObj(U,Kmap,B,O,p,lambda);
%     his = [his;fobj];
%     if(~mod(iter,100)),
%         alpha = alpha * 2;
%         beta = beta * 2;
%     end
    
end

plot(his)
ddd
% fprintf('\n');
function [fobj] = computeObj(U,O,Amap,p,Z,B,R,S,piz,pir,pis,lambda,alpha,beta)
% sum_c=1^3 lambda || r s ||_{p,1} + sum_c=1^3 <o,|z|>
%                                  + <piz, Au-b-z> + 0.5 alpha |Au-b-z|_2^2
%                                  + <pir, Dxu-r> + 0.5 beta |Dxu-r|_2^2
%                                  + <pis, Dyu-s> + 0.5 beta |Dyu-s|_2^2,

diff1 = Amap(U)-B-Z;
diff2 = difX(U)-R;
diff3 = difY(U)-S;


fobj = lambda * sum(sum((R.^p + S.^p).^(1/p))) + sum(sum(O.*abs(Z))) ...
    + mdot(piz,diff1) + 0.5*alpha*mdot(diff1,diff1)^2 ...
    + mdot(pir,diff2) + 0.5*beta *mdot(diff2,diff2)^2 ...
    + mdot(pis,diff3) + 0.5*beta *mdot(diff3,diff3)^2;


function [fobj] = computeTrueObj(U,Kmap,B,O,p,lambda)
fobj =  lambda * sum(sum((difX(U).^p + difY(U).^p).^(1/p)))  + sum(sum(O.*abs(Kmap(U)-B)));

