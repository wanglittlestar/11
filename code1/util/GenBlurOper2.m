function [P] = GenBlurOper2
% P = double(fspecial('gaussian',[7,7],5));
% P = P/sum(P(:));
 P = fspecial('average',5)
 P = P/sum(P(:));