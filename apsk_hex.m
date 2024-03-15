function y = apsk_hex(M, Es, gamma)
%M is the order of the constellation, Ec is the average constellation
%energy % n number of symbols per ring

symbols_per_ring = M/gamma;

dth = 2*pi*gamma/M; %minimum angular distance
dr = sqrt(12*Es/(4*gamma^2-1)); %minimum radial distance
phaseOff = dth/2;
y = [];
for q=1:gamma
    for p=1:symbols_per_ring
        if mod(q,2) == 0
            temp = dr/2*(2*q-1)*exp(1j*dth/2*(2*p-1))*exp(1j*phaseOff);
        else
            temp = dr/2*(2*q-1)*exp(1j*dth/2*(2*p-1));
        end
       
        y = cat(2,y,temp);
    end
end

end

