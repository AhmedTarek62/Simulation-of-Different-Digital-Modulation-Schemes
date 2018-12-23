function [seq] = quadseq(L,a,b,c,d)
seq=zeros(2,L);
for i=1:2
    for j=1:L
        temp=rand();
        if(temp<0.25)
            seq(i,j)=a;
        elseif (temp<0.5)
            seq(i,j)=b;
        elseif (temp<0.75)
            seq(i,j)=c;
        elseif (temp<1)
            seq(i,j)=d;
        end
    end
end
end

