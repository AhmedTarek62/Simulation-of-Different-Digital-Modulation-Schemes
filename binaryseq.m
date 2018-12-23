function [seq] = binaryseq(L,x,y)
seq=zeros(1,L);
for i=1:L
if(rand()>0.5)
    seq(i)=x;
else
    seq(i)=y;
end
end
end

