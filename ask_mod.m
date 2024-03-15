function symbols = ask_mod(M) %By default Es=1
%ASK_MOD Summary of this function goes here
%   Detailed explanation goes here
d = sqrt(1/((M+1)*(2*M+1)/6-1));
for i=1:M
    symbols(i) = (i-1)*d + 1j*0;
end
end

