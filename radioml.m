clear;
clc;
clear;
clc;

%Modulation orders (changes when needed)


Es = 1; %average constellation energy                           
step_snr = 5;
SNR_dB =  0:step_snr:30; %%contains 7 values
No = Es*10.^(-SNR_dB/10);
sigma_phi = sqrt(0.1); %phase noise standard deviation if the channel is GPN
kappa = 10; %Von Mises parameter should take values [1,2,3,10] k=0=>uniform distribution

num_images = 1120; %number of images per constellation frame overall snr
N=1000; %symbols per constellation frame
% ASK modulation with M=[4,8]
ask_4 = ask_mod(4);
ask_8 = ask_mod(8);

%BPSK modulation
bpsk = pskmod([0 1],2);

%QPSK modulation
qpsk = pskmod([0 1 2 3],4);

%PSK modulation with M=[8,16,32]
% psk_8 = pskmod(linspace(0,7,8),8);
% psk_16 = pskmod(linspace(0,15,16),16);
% psk_32 = pskmod(linspace(0,31,32),32);

%HQAM modulation with M=[8,16,32]
hqam_4 = hqam(4);
hqam_16 = hqam(16);
hqam_64 = hqam(64);

%QAM modulation with M=[16,32,64,128,256]
qam_16 = qammod(0:16-1,16,'UnitAveragePower',true);
qam_32 = qammod(0:32-1,32,'UnitAveragePower',true);
qam_64 = qammod(0:64-1,64,'UnitAveragePower',true);
qam_128 = qammod(0:128-1,128,'UnitAveragePower',true);
qam_256 = qammod(0:256-1,256,'UnitAveragePower',true);

%APSK modulation with M=[16,32,64,128]
apsk_16 = apskmod(linspace(0,15,16),[4 12],[0.3780 1.1339]);
apsk_32 = apskmod(linspace(0,31,32),[4 12 16],[0.25 0.75 1.25]);
apsk_64 = apskmod(linspace(0,63,64),[4 12 16 32],[0.1754 0.5262 0.8770 1.2278]);
apsk_128 = apskmod(linspace(0,127,128),[4 12 16 32 64],[0.1327 0.3982 0.6638 0.9293 1.1947]);

%By (modulation name)_order_polar we denote the modulaton in the coordinate
%system...

for i = 1:length(No)
    for k=1:num_images/length(No)

        sigma = 0.5;
        x_max = max(abs(qam_256)); %sigma and x_max define the frame size...

        % TX symbols 
        tx_ask_4 = ask_4(randi([0 3],1,N)+1); 
        tx_ask_8 = ask_8(randi([0 7],1,N)+1); 
        tx_bpsk = bpsk(randi([0 1],1,N)+1);
        tx_qpsk = qpsk(randi([0 3],1,N)+1);
%         tx_psk_8 = psk_8(randi([0 7],1,N)+1);
%         tx_psk_16 = psk_16(randi([0 15],1,N)+1);
%         tx_psk_32 = psk_32(randi([0 31],1,N)+1);
        tx_hqam_4 = hqam_4(randi([0 3],1,N)+1);
        tx_hqam_16 = hqam_16(randi([0 15],1,N)+1);
        tx_hqam_64 = hqam_64(randi([0 63],1,N)+1);
        tx_qam_16 = qam_16(randi([0 15],1,N)+1);
        tx_qam_32 = qam_32(randi([0 31],1,N)+1);
        tx_qam_64 = qam_64(randi([0 63],1,N)+1);
        tx_qam_128 = qam_128(randi([0 127],1,N)+1);
        tx_qam_256 = qam_256(randi([0 255],1,N)+1); 
        tx_apsk_16 = apsk_16(randi([0 15],1,N)+1);
        tx_apsk_32 = apsk_32(randi([0 31],1,N)+1);
        tx_apsk_64 = apsk_64(randi([0 63],1,N)+1);
        tx_apsk_128 = apsk_128(randi([0 127],1,N)+1);
        
        %Channel impairments AWGN+PN
        noise = sqrt(No(i)/2).*(randn(1,N)+1j*randn(1,N));
%        phase_noise = normrnd(0, sigma_phi , 1, N);
%         phase_noise = (rand(1,N) * 2 - 1) * pi;
%        phase_noise = vm_rand(0,kappa,[1,N]);
        phase_noise = zeros(1,N); %No phase noise
        
        %ASK sims
        rx_ask_4 = tx_ask_4.*exp(1j*phase_noise) + noise;
        rx_ask_8 = tx_ask_8.*exp(1j*phase_noise) + noise;
        [rx_ask_4_th,rx_ask_4_rho] = cart2pol(real(rx_ask_4),imag(rx_ask_4));
        [rx_ask_8_th,rx_ask_8_rho] = cart2pol(real(rx_ask_8),imag(rx_ask_8));

        rx_ask_4_polar = zeros(1,N);
        for w=1:N
            if (rx_ask_4_th(1,w) > pi)
                 rx_ask_4_polar(1,w) = rx_ask_4_rho(1,w)+1j*(rx_ask_4_th(1,w)-2*pi);
            elseif (rx_ask_4_th(1,w) < -pi)
                 rx_ask_4_polar(1,w) = rx_ask_4_rho(1,w)+1j*(rx_ask_4_th(1,w)+2*pi);
            else
                 rx_ask_4_polar(1,w) = rx_ask_4_rho(1,w)+1j*(rx_ask_4_th(1,w));
            end
        end
        rx_ask_8_polar = zeros(1,N);
        for w=1:N
            if (rx_ask_8_th(1,w) > pi)
                 rx_ask_8_polar(1,w) = rx_ask_8_rho(1,w)+1j*(rx_ask_8_th(1,w)-2*pi);
            elseif (rx_ask_8_th(1,w) < -pi)
                 rx_ask_8_polar(1,w) = rx_ask_8_rho(1,w)+1j*(rx_ask_8_th(1,w)+2*pi);
            else
                 rx_ask_8_polar(1,w) = rx_ask_8_rho(1,w)+1j*(rx_ask_8_th(1,w));
            end
        end
        
        %BPSK, QPSK sims
        rx_bpsk = tx_bpsk.*exp(1j*phase_noise) + noise;
        rx_qpsk = tx_qpsk.*exp(1j*phase_noise) + noise;
        [rx_bpsk_th,rx_bpsk_rho] = cart2pol(real(rx_bpsk),imag(rx_bpsk));
        [rx_qpsk_th,rx_qpsk_rho] = cart2pol(real(rx_qpsk),imag(rx_qpsk));

        rx_bpsk_polar = zeros(1,N);
        for w=1:N
            if (rx_bpsk_th(1,w) > pi)
                 rx_bpsk_polar(1,w) = rx_bpsk_rho(1,w)+1j*(rx_bpsk_th(1,w)-2*pi);
            elseif (rx_bpsk_th(1,w) < -pi)
                 rx_bpsk_polar(1,w) = rx_bpsk_rho(1,w)+1j*(rx_bpsk_th(1,w)+2*pi);
            else
                 rx_bpsk_polar(1,w) = rx_bpsk_rho(1,w)+1j*(rx_bpsk_th(1,w));
            end
        end
        rx_qpsk_polar = zeros(1,N);
        for w=1:N
            if (rx_qpsk_th(1,w) > pi)
                 rx_qpsk_polar(1,w) = rx_qpsk_rho(1,w)+1j*(rx_qpsk_th(1,w)-2*pi);
            elseif (rx_qpsk_th(1,w) < -pi)
                 rx_qpsk_polar(1,w) = rx_qpsk_rho(1,w)+1j*(rx_qpsk_th(1,w)+2*pi);
            else
                 rx_qpsk_polar(1,w) = rx_qpsk_rho(1,w)+1j*(rx_qpsk_th(1,w));
            end
        end
        
        %PSK sims
%         rx_hqam_4 = tx_psk_8.*exp(1j*phase_noise) + noise;
%         rx_hqam_16 = tx_psk_16.*exp(1j*phase_noise) + noise;
%         rx_hqam_64 = tx_psk_32.*exp(1j*phase_noise) + noise;
%         [rx_hqam_4_th,rx_hqam_4_rho] = cart2pol(real(rx_hqam_4),imag(rx_hqam_4));
%         [rx_hqam_16_th,rx_hqam_16_rho] = cart2pol(real(rx_hqam_16),imag(rx_hqam_16));
%         [rx_hqam_64_th,rx_hqam_64_rho] = cart2pol(real(rx_hqam_64),imag(rx_hqam_64));
%         
%         rx_hqam_4_polar = zeros(1,N);
%         for w=1:N
%             if (rx_hqam_4_th(1,w) > pi)
%                  rx_hqam_4_polar(1,w) = rx_hqam_4_rho(1,w)+1j*(rx_hqam_4_th(1,w)-2*pi);
%             elseif (rx_hqam_4_th(1,w) < -pi)
%                  rx_hqam_4_polar(1,w) = rx_hqam_4_rho(1,w)+1j*(rx_hqam_4_th(1,w)+2*pi);
%             else
%                  rx_hqam_4_polar(1,w) = rx_hqam_4_rho(1,w)+1j*(rx_hqam_4_th(1,w));
%             end
%         end
%         
%         rx_hqam_16_polar = zeros(1,N);
%         for w=1:N
%             if (rx_hqam_16_th(1,w) > pi)
%                  rx_hqam_16_polar(1,w) = rx_hqam_16_rho(1,w)+1j*(rx_hqam_16_th(1,w)-2*pi);
%             elseif (rx_hqam_16_th(1,w) < -pi)
%                  rx_hqam_16_polar(1,w) = rx_hqam_16_rho(1,w)+1j*(rx_hqam_16_th(1,w)+2*pi);
%             else
%                  rx_hqam_16_polar(1,w) = rx_hqam_16_rho(1,w)+1j*(rx_hqam_16_th(1,w));
%             end
%         end
%         
%         rx_hqam_64_polar = zeros(1,N);
%         for w=1:N
%             if (rx_hqam_64_th(1,w) > pi)
%                  rx_hqam_64_polar(1,w) = rx_hqam_64_rho(1,w)+1j*(rx_hqam_64_th(1,w)-2*pi);
%             elseif (rx_hqam_64_th(1,w) < -pi)
%                  rx_hqam_64_polar(1,w) = rx_hqam_64_rho(1,w)+1j*(rx_hqam_64_th(1,w)+2*pi);
%             else
%                  rx_hqam_64_polar(1,w) = rx_hqam_64_rho(1,w)+1j*(rx_hqam_64_th(1,w));
%             end
%         end

        %HQAM sims
        rx_hqam_4 = tx_hqam_4.*exp(1j*phase_noise) + noise;
        rx_hqam_16 = tx_hqam_16.*exp(1j*phase_noise) + noise;
        rx_hqam_64 = tx_hqam_64.*exp(1j*phase_noise) + noise;
        [rx_hqam_4_th,rx_hqam_4_rho] = cart2pol(real(rx_hqam_4),imag(rx_hqam_4));
        [rx_hqam_16_th,rx_hqam_16_rho] = cart2pol(real(rx_hqam_16),imag(rx_hqam_16));
        [rx_hqam_64_th,rx_hqam_64_rho] = cart2pol(real(rx_hqam_64),imag(rx_hqam_64));
        
        rx_hqam_4_polar = zeros(1,N);
        for w=1:N
            if (rx_hqam_4_th(1,w) > pi)
                 rx_hqam_4_polar(1,w) = rx_hqam_4_rho(1,w)+1j*(rx_hqam_4_th(1,w)-2*pi);
            elseif (rx_hqam_4_th(1,w) < -pi)
                 rx_hqam_4_polar(1,w) = rx_hqam_4_rho(1,w)+1j*(rx_hqam_4_th(1,w)+2*pi);
            else
                 rx_hqam_4_polar(1,w) = rx_hqam_4_rho(1,w)+1j*(rx_hqam_4_th(1,w));
            end
        end
        
        rx_hqam_16_polar = zeros(1,N);
        for w=1:N
            if (rx_hqam_16_th(1,w) > pi)
                 rx_hqam_16_polar(1,w) = rx_hqam_16_rho(1,w)+1j*(rx_hqam_16_th(1,w)-2*pi);
            elseif (rx_hqam_16_th(1,w) < -pi)
                 rx_hqam_16_polar(1,w) = rx_hqam_16_rho(1,w)+1j*(rx_hqam_16_th(1,w)+2*pi);
            else
                 rx_hqam_16_polar(1,w) = rx_hqam_16_rho(1,w)+1j*(rx_hqam_16_th(1,w));
            end
        end
        
        rx_hqam_64_polar = zeros(1,N);
        for w=1:N
            if (rx_hqam_64_th(1,w) > pi)
                 rx_hqam_64_polar(1,w) = rx_hqam_64_rho(1,w)+1j*(rx_hqam_64_th(1,w)-2*pi);
            elseif (rx_hqam_64_th(1,w) < -pi)
                 rx_hqam_64_polar(1,w) = rx_hqam_64_rho(1,w)+1j*(rx_hqam_64_th(1,w)+2*pi);
            else
                 rx_hqam_64_polar(1,w) = rx_hqam_64_rho(1,w)+1j*(rx_hqam_64_th(1,w));
            end
        end
        
        %QAM sims
        rx_qam_16 = tx_qam_16.*exp(1j*phase_noise) + noise;
        rx_qam_32 = tx_qam_32.*exp(1j*phase_noise) + noise;
        rx_qam_64= tx_qam_64.*exp(1j*phase_noise) + noise;
        rx_qam_128 = tx_qam_128.*exp(1j*phase_noise) + noise;
        rx_qam_256 = tx_qam_256.*exp(1j*phase_noise) + noise;
        [rx_qam_16_th,rx_qam_16_rho] = cart2pol(real(rx_qam_16),imag(rx_qam_16));
        [rx_qam_32_th,rx_qam_32_rho] = cart2pol(real(rx_qam_32),imag(rx_qam_32));
        [rx_qam_64_th,rx_qam_64_rho] = cart2pol(real(rx_qam_64),imag(rx_qam_64));
        [rx_qam_128_th,rx_qam_128_rho] = cart2pol(real(rx_qam_128),imag(rx_qam_128));
        [rx_qam_256_th,rx_qam_256_rho] = cart2pol(real(rx_qam_256),imag(rx_qam_256));
        
        rx_qam_16_polar = zeros(1,N);
        for w=1:N
            if (rx_qam_16_th(1,w) > pi)
                 rx_qam_16_polar(1,w) = rx_qam_16_rho(1,w)+1j*(rx_qam_16_th(1,w)-2*pi);
            elseif (rx_qam_16_th(1,w) < -pi)
                 rx_qam_16_polar(1,w) = rx_qam_16_rho(1,w)+1j*(rx_qam_16_th(1,w)+2*pi);
            else
                 rx_qam_16_polar(1,w) = rx_qam_16_rho(1,w)+1j*(rx_qam_16_th(1,w));
            end
        end
        
        rx_qam_32_polar = zeros(1,N);
        for w=1:N
            if (rx_qam_32_th(1,w) > pi)
                 rx_qam_32_polar(1,w) = rx_qam_32_rho(1,w)+1j*(rx_qam_32_th(1,w)-2*pi);
            elseif (rx_qam_32_th(1,w) < -pi)
                 rx_qam_32_polar(1,w) = rx_qam_32_rho(1,w)+1j*(rx_qam_32_th(1,w)+2*pi);
            else
                 rx_qam_32_polar(1,w) = rx_qam_32_rho(1,w)+1j*(rx_qam_32_th(1,w));
            end
        end
        
        rx_qam_64_polar = zeros(1,N);
        for w=1:N
            if (rx_qam_64_th(1,w) > pi)
                 rx_qam_64_polar(1,w) = rx_qam_64_rho(1,w)+1j*(rx_qam_64_th(1,w)-2*pi);
            elseif (rx_qam_64_th(1,w) < -pi)
                 rx_qam_64_polar(1,w) = rx_qam_64_rho(1,w)+1j*(rx_qam_64_th(1,w)+2*pi);
            else
                 rx_qam_64_polar(1,w) = rx_qam_64_rho(1,w)+1j*(rx_qam_64_th(1,w));
            end
        end
        
        rx_qam_128_polar = zeros(1,N);
        for w=1:N
            if (rx_qam_128_th(1,w) > pi)
                 rx_qam_128_polar(1,w) = rx_qam_128_rho(1,w)+1j*(rx_qam_128_th(1,w)-2*pi);
            elseif (rx_qam_128_th(1,w) < -pi)
                 rx_qam_128_polar(1,w) = rx_qam_128_rho(1,w)+1j*(rx_qam_128_th(1,w)+2*pi);
            else
                 rx_qam_128_polar(1,w) = rx_qam_128_rho(1,w)+1j*(rx_qam_128_th(1,w));
            end
        end
        
        rx_qam_256_polar = zeros(1,N);
        for w=1:N
            if (rx_qam_256_th(1,w) > pi)
                 rx_qam_256_polar(1,w) = rx_qam_256_rho(1,w)+1j*(rx_qam_256_th(1,w)-2*pi);
            elseif (rx_qam_256_th(1,w) < -pi)
                 rx_qam_256_polar(1,w) = rx_qam_256_rho(1,w)+1j*(rx_qam_256_th(1,w)+2*pi);
            else
                 rx_qam_256_polar(1,w) = rx_qam_256_rho(1,w)+1j*(rx_qam_256_th(1,w));
            end
        end
        
        %APSK sims
        rx_apsk_16 = tx_apsk_16.*exp(1j*phase_noise) + noise;
        rx_apsk_32 = tx_apsk_32.*exp(1j*phase_noise) + noise;
        rx_apsk_64 = tx_apsk_64.*exp(1j*phase_noise) + noise;
        rx_apsk_128 = tx_apsk_128.*exp(1j*phase_noise) + noise;
        [rx_apsk_16_th,rx_apsk_16_rho] = cart2pol(real(rx_apsk_16),imag(rx_apsk_16));
        [rx_apsk_32_th,rx_apsk_32_rho] = cart2pol(real(rx_apsk_32),imag(rx_apsk_32));
        [rx_apsk_64_th,rx_apsk_64_rho] = cart2pol(real(rx_apsk_64),imag(rx_apsk_64));
        [rx_apsk_128_th,rx_apsk_128_rho] = cart2pol(real(rx_apsk_128),imag(rx_apsk_128));
       
        
        rx_apsk_16_polar = zeros(1,N);
        for w=1:N
            if (rx_apsk_16_th(1,w) > pi)
                 rx_apsk_16_polar(1,w) = rx_apsk_16_rho(1,w)+1j*(rx_apsk_16_th(1,w)-2*pi);
            elseif (rx_apsk_16_th(1,w) < -pi)
                 rx_apsk_16_polar(1,w) = rx_apsk_16_rho(1,w)+1j*(rx_apsk_16_th(1,w)+2*pi);
            else
                 rx_apsk_16_polar(1,w) = rx_apsk_16_rho(1,w)+1j*(rx_apsk_16_th(1,w));
            end
        end
        
        rx_apsk_32_polar = zeros(1,N);
        for w=1:N
            if (rx_apsk_32_th(1,w) > pi)
                 rx_apsk_32_polar(1,w) = rx_apsk_32_rho(1,w)+1j*(rx_apsk_32_th(1,w)-2*pi);
            elseif (rx_apsk_32_th(1,w) < -pi)
                 rx_apsk_32_polar(1,w) = rx_apsk_32_rho(1,w)+1j*(rx_apsk_32_th(1,w)+2*pi);
            else
                 rx_apsk_32_polar(1,w) = rx_apsk_32_rho(1,w)+1j*(rx_apsk_32_th(1,w));
            end
        end
        
        rx_apsk_64_polar = zeros(1,N);
        for w=1:N
            if (rx_apsk_64_th(1,w) > pi)
                 rx_apsk_64_polar(1,w) = rx_apsk_64_rho(1,w)+1j*(rx_apsk_64_th(1,w)-2*pi);
            elseif (rx_apsk_64_th(1,w) < -pi)
                 rx_apsk_64_polar(1,w) = rx_apsk_64_rho(1,w)+1j*(rx_apsk_64_th(1,w)+2*pi);
            else
                 rx_apsk_64_polar(1,w) = rx_apsk_64_rho(1,w)+1j*(rx_apsk_64_th(1,w));
            end
        end
        
        rx_apsk_128_polar = zeros(1,N);
        for w=1:N
            if (rx_apsk_128_th(1,w) > pi)
                 rx_apsk_128_polar(1,w) = rx_apsk_128_rho(1,w)+1j*(rx_apsk_128_th(1,w)-2*pi);
            elseif (rx_apsk_128_th(1,w) < -pi)
                 rx_apsk_128_polar(1,w) = rx_apsk_128_rho(1,w)+1j*(rx_apsk_128_th(1,w)+2*pi);
            else
                 rx_apsk_128_polar(1,w) = rx_apsk_128_rho(1,w)+1j*(rx_apsk_128_th(1,w));
            end
        end

        % Set axis limits 
        x_axis_limit_iq = [-sqrt(2)*(x_max+sigma), sqrt(2)*(x_max+sigma)];
        y_axis_limit_iq = [-sqrt(2)*(x_max+sigma), sqrt(2)*(x_max+sigma)];
        x_axis_limit_polar = [0, sqrt(2)*(x_max+sigma)];
        y_axis_limit_polar = [-pi , pi];
        
        %ASK images IQ
        image_ask_4 = plot(real(rx_ask_4) , imag(rx_ask_4) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_iq);
        ylim(y_axis_limit_iq);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_ask_4{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_ask_4{((i-1)*num_images/length(No))+k} = imresize(images_ask_4{((i-1)*num_images/length(No))+k},2/3); 
        
        image_ask_8 = plot(real(rx_ask_8) , imag(rx_ask_8) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_iq);
        ylim(y_axis_limit_iq);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_ask_8{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_ask_8{((i-1)*num_images/length(No))+k} = imresize(images_ask_8{((i-1)*num_images/length(No))+k},2/3);
    
        %ASK images POLAR
        image_ask_4_polar = plot(real(rx_ask_4_polar) , imag(rx_ask_4_polar) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_polar);
        ylim(y_axis_limit_polar);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_ask_4_polar{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_ask_4_polar{((i-1)*num_images/length(No))+k} = imresize(images_ask_4_polar{((i-1)*num_images/length(No))+k},2/3);

        image_ask_8_polar = plot(real(rx_ask_8_polar) , imag(rx_ask_8_polar) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_polar);
        ylim(y_axis_limit_polar);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_ask_8_polar{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_ask_8_polar{((i-1)*num_images/length(No))+k} = imresize(images_ask_8_polar{((i-1)*num_images/length(No))+k},2/3);
        
        %B-QPSK images IQ
        image_bpsk = plot(real(rx_bpsk) , imag(rx_bpsk) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_iq);
        ylim(y_axis_limit_iq);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_bpsk{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_bpsk{((i-1)*num_images/length(No))+k} = imresize(images_bpsk{((i-1)*num_images/length(No))+k},2/3);

        image_qpsk = plot(real(rx_qpsk) , imag(rx_qpsk) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_iq);
        ylim(y_axis_limit_iq);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_qpsk{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_qpsk{((i-1)*num_images/length(No))+k} = imresize(images_qpsk{((i-1)*num_images/length(No))+k},2/3);
    
        %ASK images POLAR
        image_bpsk_polar = plot(real(rx_bpsk_polar) , imag(rx_bpsk_polar) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_polar);
        ylim(y_axis_limit_polar);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_bpsk_polar{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_bpsk_polar{((i-1)*num_images/length(No))+k} = imresize(images_bpsk_polar{((i-1)*num_images/length(No))+k},2/3);
        
        image_qpsk_polar = plot(real(rx_qpsk_polar) , imag(rx_qpsk_polar) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_polar);
        ylim(y_axis_limit_polar);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_qpsk_polar{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_qpsk_polar{((i-1)*num_images/length(No))+k} = imresize(images_qpsk_polar{((i-1)*num_images/length(No))+k},2/3);

        %HQAM images IQ
        image_hqam_4 = plot(real(rx_hqam_4) , imag(rx_hqam_4) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_iq);
        ylim(y_axis_limit_iq);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_hqam_4{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_psk_8{((i-1)*num_images/length(No))+k} = imresize(images_psk_8{((i-1)*num_images/length(No))+k},2/3);

        image_hqam_16 = plot(real(rx_hqam_16) , imag(rx_hqam_16) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_iq);
        ylim(y_axis_limit_iq);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_hqam_16{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_psk_16{((i-1)*num_images/length(No))+k} = imresize(images_psk_16{((i-1)*num_images/length(No))+k},2/3);

        image_hqam_64 = plot(real(rx_hqam_64) , imag(rx_hqam_64) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_iq);
        ylim(y_axis_limit_iq);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_hqam_64{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_psk_32{((i-1)*num_images/length(No))+k} = imresize(images_psk_32{((i-1)*num_images/length(No))+k},2/3);
        
        %HQAM images polar
        image_hqam_4_polar = plot(real(rx_hqam_4_polar) , imag(rx_hqam_4_polar) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_polar);
        ylim(y_axis_limit_polar);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_hqam_4_polar{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_psk_8_polar{((i-1)*num_images/length(No))+k} = imresize(images_psk_8_polar{((i-1)*num_images/length(No))+k},2/3);

        image_hqam_16_polar = plot(real(rx_hqam_16_polar) , imag(rx_hqam_16_polar) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_polar);
        ylim(y_axis_limit_polar);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_hqam_16_polar{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_psk_16_polar{((i-1)*num_images/length(No))+k} = imresize(images_psk_16_polar{((i-1)*num_images/length(No))+k},2/3);

        image_hqam_64_polar = plot(real(rx_hqam_64_polar) , imag(rx_hqam_64_polar) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_polar);
        ylim(y_axis_limit_polar);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_hqam_64_polar{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_psk_32_polar{((i-1)*num_images/length(No))+k} = imresize(images_psk_32_polar{((i-1)*num_images/length(No))+k},2/3);

        %QAM images IQ
        image_qam_16 = plot(real(rx_qam_16) , imag(rx_qam_16) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_iq);
        ylim(y_axis_limit_iq);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_qam_16{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_qam_16{((i-1)*num_images/length(No))+k} = imresize(images_qam_16{((i-1)*num_images/length(No))+k},2/3);

        image_qam_32 = plot(real(rx_qam_32) , imag(rx_qam_32) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_iq);
        ylim(y_axis_limit_iq);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_qam_32{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_qam_32{((i-1)*num_images/length(No))+k} = imresize(images_qam_32{((i-1)*num_images/length(No))+k},2/3);

        image_qam_64 = plot(real(rx_qam_64) , imag(rx_qam_64) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_iq);
        ylim(y_axis_limit_iq);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_qam_64{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_qam_64{((i-1)*num_images/length(No))+k} = imresize(images_qam_64{((i-1)*num_images/length(No))+k},2/3);

        image_qam_128 = plot(real(rx_qam_128) , imag(rx_qam_128) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_iq);
        ylim(y_axis_limit_iq);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_qam_128{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_qam_128{((i-1)*num_images/length(No))+k} = imresize(images_qam_128{((i-1)*num_images/length(No))+k},2/3);

        image_qam_256 = plot(real(rx_qam_256) , imag(rx_qam_256) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_iq);
        ylim(y_axis_limit_iq);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_qam_256{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_qam_256{((i-1)*num_images/length(No))+k} = imresize(images_qam_256{((i-1)*num_images/length(No))+k},2/3);

        %QAM images POLAR
        image_qam_16_polar = plot(real(rx_qam_16_polar) , imag(rx_qam_16_polar) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_polar);
        ylim(y_axis_limit_polar);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_qam_16_polar{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_qam_16_polar{((i-1)*num_images/length(No))+k} = imresize(images_qam_16_polar{((i-1)*num_images/length(No))+k},2/3);

        image_qam_32_polar = plot(real(rx_qam_32_polar) , imag(rx_qam_32_polar) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_polar);
        ylim(y_axis_limit_polar);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_qam_32_polar{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_qam_32_polar{((i-1)*num_images/length(No))+k} = imresize(images_qam_32_polar{((i-1)*num_images/length(No))+k},2/3);

        image_qam_64_polar = plot(real(rx_qam_64_polar) , imag(rx_qam_64_polar) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_polar);
        ylim(y_axis_limit_polar);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_qam_64_polar{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_qam_64_polar{((i-1)*num_images/length(No))+k} = imresize(images_qam_64_polar{((i-1)*num_images/length(No))+k},2/3);

        image_qam_128_polar = plot(real(rx_qam_128_polar) , imag(rx_qam_128_polar) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_polar);
        ylim(y_axis_limit_polar);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_qam_128_polar{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_qam_128_polar{((i-1)*num_images/length(No))+k} = imresize(images_qam_128_polar{((i-1)*num_images/length(No))+k},2/3);

        image_qam_256_polar = plot(real(rx_qam_256_polar) , imag(rx_qam_256_polar) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_polar);
        ylim(y_axis_limit_polar);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_qam_256_polar{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_qam_256_polar{((i-1)*num_images/length(No))+k} = imresize(images_qam_256_polar{((i-1)*num_images/length(No))+k},2/3);

        %APSK images IQ
        image_apsk_16 = plot(real(rx_apsk_16) , imag(rx_apsk_16) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_iq);
        ylim(y_axis_limit_iq);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_apsk_16{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_apsk_16{((i-1)*num_images/length(No))+k} = imresize(images_apsk_16{((i-1)*num_images/length(No))+k},2/3);

        image_apsk_32 = plot(real(rx_apsk_32) , imag(rx_apsk_32) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_iq);
        ylim(y_axis_limit_iq);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_apsk_32{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_apsk_32{((i-1)*num_images/length(No))+k} = imresize(images_apsk_32{((i-1)*num_images/length(No))+k},2/3);

        image_apsk_64 = plot(real(rx_apsk_64) , imag(rx_apsk_64) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_iq);
        ylim(y_axis_limit_iq);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_apsk_64{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_apsk_64{((i-1)*num_images/length(No))+k} = imresize(images_apsk_64{((i-1)*num_images/length(No))+k},2/3);

        image_apsk_128 = plot(real(rx_apsk_128) , imag(rx_apsk_128) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_iq);
        ylim(y_axis_limit_iq);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_apsk_128{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_apsk_128{((i-1)*num_images/length(No))+k} = imresize(images_apsk_128{((i-1)*num_images/length(No))+k},2/3);

        %APSK images POLAR
        image_apsk_16_polar = plot(real(rx_apsk_16_polar) , imag(rx_apsk_16_polar) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_polar);
        ylim(y_axis_limit_polar);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_apsk_16_polar{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_apsk_16_polar{((i-1)*num_images/length(No))+k} = imresize(images_apsk_16_polar{((i-1)*num_images/length(No))+k},2/3);

        image_apsk_32_polar = plot(real(rx_apsk_32_polar) , imag(rx_apsk_32_polar) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_polar);
        ylim(y_axis_limit_polar);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_apsk_32_polar{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_apsk_32_polar{((i-1)*num_images/length(No))+k} = imresize(images_apsk_32_polar{((i-1)*num_images/length(No))+k},2/3);
        
        image_apsk_64_polar = plot(real(rx_apsk_64_polar) , imag(rx_apsk_64_polar) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_polar);
        ylim(y_axis_limit_polar);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_apsk_64_polar{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_apsk_64_polar{((i-1)*num_images/length(No))+k} = imresize(images_apsk_64_polar{((i-1)*num_images/length(No))+k},2/3);
        
        image_apsk_128_polar = plot(real(rx_apsk_128_polar) , imag(rx_apsk_128_polar) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_polar);
        ylim(y_axis_limit_polar);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_apsk_128_polar{((i-1)*num_images/length(No))+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        %images_apsk_128_polar{((i-1)*num_images/length(No))+k} = imresize(images_apsk_128_polar{((i-1)*num_images/length(No))+k},2/3);
    end
end



%Create label matrix
labels_before = [];
for m = 1:16
    labels_before = [labels_before; m*ones(num_images, 1)];
end


num_classes = 16;
labels_onehot_before = zeros(length(labels_before), num_classes);

for i = 1:length(labels_before)
    labels_onehot_before(i, labels_before(i)) = 1;
end


labels_onehot_total_before = labels_onehot_before;

%Create label matrix
labels_after = [];
for m = 1:16
    labels_after = [labels_after; m*ones(num_images, 1)];
end

labels_onehot_after = zeros(length(labels_after), num_classes);

for i = 1:length(labels_after)
    labels_onehot_after(i, labels_after(i)) = 1;
end

labels_onehot_total_after = labels_onehot_after;


%Convert the cell array to a regular array
Images_ask_4 = cell2mat(images_ask_4);
Images_ask_8 = cell2mat(images_ask_8);

Images_bpsk = cell2mat(images_bpsk);
Images_qpsk = cell2mat(images_qpsk);

% Images_psk_8 = cell2mat(images_psk_8);
% Images_psk_16 = cell2mat(images_psk_16);
% Images_psk_32 = cell2mat(images_psk_32);

Images_hqam_4 = cell2mat(images_hqam_4);
Images_hqam_16 = cell2mat(images_hqam_16);
Images_hqam_64 = cell2mat(images_hqam_64);

Images_qam_16 = cell2mat(images_qam_16);
Images_qam_32 = cell2mat(images_qam_32);
Images_qam_64 = cell2mat(images_qam_64);
Images_qam_128 = cell2mat(images_qam_128);
Images_qam_256 = cell2mat(images_qam_256);

Images_apsk_16 = cell2mat(images_apsk_16);
Images_apsk_32 = cell2mat(images_apsk_32);
Images_apsk_64 = cell2mat(images_apsk_64);
Images_apsk_128 = cell2mat(images_apsk_128);

% Convert polar forms
Images_ask_4_polar = cell2mat(images_ask_4_polar);
Images_ask_8_polar = cell2mat(images_ask_8_polar);

Images_bpsk_polar = cell2mat(images_bpsk_polar);
Images_qpsk_polar = cell2mat(images_qpsk_polar);

% Images_psk_8_polar = cell2mat(images_psk_8_polar);
% Images_psk_16_polar = cell2mat(images_psk_16_polar);
% Images_psk_32_polar = cell2mat(images_psk_32_polar);

Images_hqam_4_polar = cell2mat(images_hqam_4_polar);
Images_hqam_16_polar = cell2mat(images_hqam_16_polar);
Images_hqam_64_polar = cell2mat(images_hqam_64_polar);

Images_qam_16_polar = cell2mat(images_qam_16_polar);
Images_qam_32_polar = cell2mat(images_qam_32_polar);
Images_qam_64_polar = cell2mat(images_qam_64_polar);
Images_qam_128_polar = cell2mat(images_qam_128_polar);
Images_qam_256_polar = cell2mat(images_qam_256_polar);

Images_apsk_16_polar = cell2mat(images_apsk_16_polar);
Images_apsk_32_polar = cell2mat(images_apsk_32_polar);
Images_apsk_64_polar = cell2mat(images_apsk_64_polar);
Images_apsk_128_polar = cell2mat(images_apsk_128_polar);

% Concatenated Images
Images_total_before = [Images_ask_4'; Images_ask_8'; Images_bpsk'; Images_qpsk'; Images_hqam_4'; Images_hqam_16'; ...
    Images_hqam_64'; Images_qam_16'; Images_qam_32'; Images_qam_64'; Images_qam_128'; Images_qam_256'; ...
    Images_apsk_16'; Images_apsk_32'; Images_apsk_64'; Images_apsk_128';];
Images_total_after = [Images_ask_4_polar'; Images_ask_8_polar'; Images_bpsk_polar'; Images_qpsk_polar'; ...
    Images_hqam_4_polar'; Images_hqam_16_polar'; Images_hqam_64_polar'; Images_qam_16_polar'; Images_qam_32_polar'; ...
    Images_qam_64_polar'; Images_qam_128_polar'; Images_qam_256_polar'; Images_apsk_16_polar'; Images_apsk_32_polar'; ...
    Images_apsk_64_polar'; Images_apsk_128_polar';];

Image_order = zeros(1, num_images);
idx = 1;

for snr=min(SNR_dB):step_snr:max(SNR_dB)
    for i = 1:(num_images/length(No))
        Image_order(idx) = snr;
        idx = idx + 1;
    end
end

% Image_order = transpose(image_order);
% 
% image_order_v2 = zeros(1, num_images*length_M);
% idx = 1;
% for m=1:length_M
%     for snr=min(SNR_dB):step_snr:max(SNR_dB)
%         for i = 1:(num_images/length(No))
%             image_order_v2(idx) = snr;
%             idx = idx + 1;
%         end
%     end
% end
% Image_order_v2 = transpose(image_order_v2);
% 
% image_order_q = zeros(1, num_images*length_M);
% idx = 1;
% for m=1:length_M
%     for snr=min(SNR_dB):step_snr:max(SNR_dB)
%         for i = 1:(num_images/length(No))
%             image_order_q(idx) = snr;
%             idx = idx + 1;
%         end
%     end
% end
% Image_order_q = transpose(image_order_q);
% 
% image_order_q_v2 = zeros(1, num_images*length_M);
% idx = 1;
% for m=1:length_M
%     for snr=min(SNR_dB):step_snr:max(SNR_dB)
%         for i = 1:(num_images/length(No))
%             image_order_q_v2(idx) = snr;
%             idx = idx + 1;
%         end
%     end
% end
% Image_order_q_v2 = transpose(image_order_q_v2);
% 
% 
% image_order_qam = zeros(1, num_images*length_M);
% idx = 1;
% for m=1:length_M
%     for snr=min(SNR_dB):step_snr:max(SNR_dB)
%         for i = 1:(num_images/length(No))
%             image_order_qam(idx) = snr;
%             idx = idx + 1;
%         end
%     end
% end
% Image_order_qam = transpose(image_order_qam);
% 
% image_order_qam_v2 = zeros(1, num_images*length_M);
% idx = 1;
% for m=1:length_M
%     for snr=min(SNR_dB):step_snr:max(SNR_dB)
%         for i = 1:(num_images/length(No))
%             image_order_qam_v2(idx) = snr;
%             idx = idx + 1;
%         end
%     end
% end
% Image_order_qam_v2 = transpose(image_order_qam_v2);
% 
% 
% image_order_hqam = zeros(1, num_images*length_M);
% idx = 1;
% for m=1:length_M
%     for snr=min(SNR_dB):step_snr:max(SNR_dB)
%         for i = 1:(num_images/length(No))
%             image_order_hqam(idx) = snr;
%             idx = idx + 1;
%         end
%     end
% end
% Image_order_hqam = transpose(image_order_hqam);
% 
% image_order_hqam_v2 = zeros(1, num_images*length_M);
% idx = 1;
% for m=1:length_M
%     for snr=min(SNR_dB):step_snr:max(SNR_dB)
%         for i = 1:(num_images/length(No))
%             image_order_hqam_v2(idx) = snr;
%             idx = idx + 1;
%         end
%     end
% end
% Image_order_hqam_v2 = transpose(image_order_hqam_v2);
Image_order_total_before = [];
for ll = 1:16
    Image_order_total_before = horzcat(Image_order_total_before, Image_order);
end
Image_order_total_before = Image_order_total_before';

Image_order_total_after = [];
for ll = 1:16
    Image_order_total_after = horzcat(Image_order_total_after, Image_order);
end
Image_order_total_after = Image_order_total_after';
% 
% Image_order_total_before = [Image_order ; Image_order_q;  Image_order_qam ; Image_order_hqam ;];
% Image_order_total_after = [Image_order_v2 ; Image_order_q_v2; Image_order_qam_v2;  Image_order_hqam_v2;];
%Modify the following command to save the updated  datasets 
%save('route' ,'Images_total_before' ,'Images_total_after','labels_onehot_total_before' , 'labels_onehot_total_after','Image_order_total_before','Image_order_total_after','-v7.3')