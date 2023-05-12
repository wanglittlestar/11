
function R = divX(Px)
R = Px-Px(:,[1 1:end-1],:);
R(:,1,:)   = Px(:,1,:);   
R(:,end,:) = -Px(:,end-1,:);
