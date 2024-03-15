function y = polar_apsk(M,Es, gamma)
symbols_per_ring = M/gamma;

dth = 2*pi*gamma/M;             %minimum angular distance
dr = sqrt(12*Es/(4*gamma^2-1)); %minimum radial distance
y = [];
for p=1:symbols_per_ring
    for q=1:gamma
        temp = dr/2*(2*q-1)*exp(1j*dth/2*(2*p-1));
        y = cat(2,y,temp);
    end
end
end

