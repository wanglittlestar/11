function [R,S] = threadholding_RS(A,B,lambda,beta,p)
% Solve the following OP:
% min_{R,S} lambda sum(sum((R.^p + S.^p).^(1/p))) + 0.5 beta ||R + A||_F^2 + 0.5 beta ||S + B||_F^2

if(p==1)
     R = - sign(A).*max(0,abs(A)-(lambda/beta));
     S = - sign(B).*max(0,abs(B)-(lambda/beta));
elseif(p==2)
     normRS = sqrt(A.^2 + B.^2);
     normRS(normRS<1e-12)=1;
     Weight = max(0,normRS-(lambda/beta))./ normRS;
     Weight(normRS==0)=1;
     R = -A.*Weight;
     S = -B.*Weight; 
end

