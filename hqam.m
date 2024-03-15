function symbols = hqam(M)
%HQAM Summary of this function goes here
%   Detailed explanation goes here
side_length = sqrt(M);
d = sqrt(12/(7*M-4));
dmin = 2*d;
symbols = zeros(1,M);
idx=0;
for i=1:side_length
    for j=1:side_length
        idx=idx+1;
        x = dmin*i+dmin*0.5*(mod(j,2)) ;
        y = dmin*sqrt(3)*j/2;
        symbols(idx) = x + 1j*y;
    end
    
end
maxx = max(real(symbols));
minx = min(real(symbols));
maxy = max(imag(symbols));
miny = min(imag(symbols));
mid_x = (maxx+minx)/2;
mid_y = (maxy+miny)/2;
for i=1:M
    symbols(i) = real(symbols(i))-mid_x + 1j*(imag(symbols(i))-mid_y);
end
end

