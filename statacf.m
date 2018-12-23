function [acf] = statacf(x)
[r,c] = size(x);
acf=zeros(1,2*c-1);
for i=1:r
    temp=conv(x(i,:),fliplr(x(i,:)));
    acf= acf+temp;
end
acf=acf/r;
end

