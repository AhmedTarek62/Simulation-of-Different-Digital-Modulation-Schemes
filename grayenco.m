function [code] = grayenco(x,a,b,c,d)
% a,b,c & d are symbol amplitudes in an ascending order
code=strings(1,length(x));
for i=1:length(x)
    if(x(1,i)==a && x(2,i)==c)
        code(i)='1010';
    elseif(x(1,i)==b && x(2,i)==c)
        code(i)='1000';
    elseif(x(1,i)==c && x(2,i)==c)
        code(i)='1100';
     elseif(x(1,i)==d && x(2,i)==c)
        code(i)='1101';
    elseif(x(1,i)==a && x(2,i)==b)
        code(i)='0001';
    elseif(x(1,i)==b && x(2,i)==b)
        code(i)='0000';
    elseif(x(1,i)==c && x(2,i)==b)
        code(i)='0100';
    elseif(x(1,i)==d && x(2,i)==b)
        code(i)='0110';
    elseif(x(1,i)==a && x(2,i)==a)
        code(i)='0011';
    elseif(x(1,i)==b && x(2,i)==a)
        code(i)='0010';
    elseif(x(1,i)==c && x(2,i)==a)
        code(i)='0101';
    elseif(x(1,i)==d && x(2,i)==a)
        code(i)='0111';
    elseif(x(1,i)==a && x(2,i)==d)
        code(i)='1011';
    elseif(x(1,i)==b && x(2,i)==d)
        code(i)='1001';
    elseif(x(1,i)==c && x(2,i)==d)
        code(i)='1110';
    elseif(x(1,i)==d && x(2,i)==d)
        code(i)='1111';
    end
end
end

