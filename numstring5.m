function [s] = numstring5(n)

if n<10
    s=strcat('0000',int2str(n));
elseif n<100
    s=strcat('000',int2str(n));
elseif n<1000
    s=strcat('00',int2str(n));
elseif n<10000
    s=strcat('0',int2str(n));
else
    s=int2str(n);
end

return
