function flag = checkdec(a)
last=inf;flag=1;
for i=1:length(a),
    if(a(i)>last && abs(a(i)-last)/abs(a(i)) >1e-6)
        fprintf('a: %e, b: %e, diff: %e\n',a(i),last,a(i)-last);
        flag=0;
        error('NonDec');
        return;
    end
    last=a(i);
end