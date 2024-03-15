clear;
clc;

%Modulation orders (changes when needed)
M = [16,32,64];

Es = 1;                             
step_snr = 5;

gamma1 = 8;                    
gamma2 = 16;

length_M = length(M);

num_images = 105;           

SNR_dB =  0:step_snr:30;
No = Es*10.^(-SNR_dB/10);

%Case 2: Channel with Phase noise and AWGN
mi_n = 0;
mi_phi = 0;

sigma_phi = 1;

N=1000; %symbols

images_length = length(No)*length_M*num_images/length(No);

images = cell(images_length, 1);
images_v2 = cell(images_length, 1);

images_q = cell(images_length, 1);
images_q_v2 = cell(images_length, 1);

images_qam = cell(images_length, 1);
images_qam_v2 = cell(images_length, 1);

images_hqam = cell(images_length, 1);
images_hqam_v2 = cell(images_length, 1);

constellations1 = [];  % Array to store constellation values for M
constellations2 = [];  % Array to store constellation values for M
constellations_qam = []; % Array to store constellation values for M
constellations_hqam = []; % Array to store constellation values for M

% Loop over the range of SNR values and modulation orders
for m = 1:length_M
  for i = 1:length(No)
    for k=1:num_images/length(No)
        
        sigma_n = sqrt(No(i));
        
        sigma = 0.1;

        gamma_squared = sigma_phi^2/sigma_n^2 + 1/Es;


        %M-APSK 
        
        
        constellation = polar_apsk(M(m),Es,gamma1); %apsk(M(m),gamma1)

        constellations1 = [constellations1 constellation];
        constellations1_unique = unique(constellations1, 'rows');

        x1 = abs(constellations1_unique);

        % Find the maximum of each row
        max_per_row_1 = max(x1, [], 2);
        
        % Find the maximum among the maximum values of each row
        x1_max = max(max_per_row_1);
        
        Constellation = polar_apsk(M(m),Es,gamma2); %apsk(M(m),gamma2)

        constellations2 = [constellations2 Constellation];

        constellations2_unique = unique(constellations2, 'rows');

        x2 = abs(constellations2_unique);
       
        % Find the maximum of each row
        max_per_row_2 = max(x2, [], 2);
        
        % Find the maximum among the maximum values of each row
        x2_max = max(max_per_row_2);
       
        Constellation_qam = qammod(0:M(m)-1,M(m),'UnitAveragePower',true); %qam(M(m))
        
        x3_max = 1.07;
        
        Constellation_hqam = hqam(M(m)); %hqam(M(m))
        
        x4_max =  1.6348;
        
        
        x_max = max([x1_max,x2_max,x3_max,x4_max]);
        % Generate random data symbols
        data_symbols = randi([0 M(m)-1],1,N); 
        Data_symbols = randi([0 M(m)-1],1,N);
        
        data_symbols_qam = randi([0 M(m)-1],1,N);
                
        data_symbols_hqam = randi([0 M(m)-1],1,N);

        
        % Map the data symbols to the M-APSK constellation points
        transmitted_symbols = constellation(data_symbols+1); %TX #N symbols from apsk(M(m),gamma1)
        Transmitted_symbols = Constellation(Data_symbols+1); %TX #N symbols from apsk(M(m),gamma2)

        transmitted_symbols_qam = Constellation_qam(data_symbols_qam+1);  %TX #N symbols from qam(M(m),gamma1)
        
        transmitted_symbols_hqam = Constellation_hqam(data_symbols_hqam+1); %TX #N symbols from qam(M(m),gamma1)
        
        
        noise = sqrt(No(i)/2).*(randn(1,N)+1j*randn(1,N));
        %phase_noise = normrnd(mi_n , sigma_phi , 1, N);
        phase_noise = (rand(1,N) * 2 - 1) * pi;
        
        %APSK(M(m),gamma1)
        s_x = real(transmitted_symbols);
        s_y = imag(transmitted_symbols);
                   
        [s_theta,s_rho] = cart2pol(s_x,s_y);       
        

        const_symbols = s_rho+1j*s_theta;

        % Generate complex Gaussian noise with zero mean and unit variance
        
        received_symbols = transmitted_symbols.*exp(1j*phase_noise) + noise;
        [r_theta,r_rho] = cart2pol(real(received_symbols),imag(received_symbols));
        

        last_symbols = zeros(1,N);
        for w=1:N
            if (r_theta(1,w) > pi)
                 last_symbols(1,w) = r_rho(1,w)+1j*(r_theta(1,w)-2*pi);
            elseif (r_theta(1,w) < -pi)
                 last_symbols(1,w) = r_rho(1,w)+1j*(r_theta(1,w)+2*pi);
            else
                 last_symbols(1,w) = r_rho(1,w)+1j*(r_theta(1,w));
            end
        end
       
        %APSK(M(m),gamma2)
        S_x = real(Transmitted_symbols);
        S_y = imag(Transmitted_symbols);

        [S_theta,S_rho] = cart2pol(S_x,S_y);
      

        Const_symbols = S_rho+1j*S_theta;
        
        Received_symbols = Transmitted_symbols.*exp(1j*phase_noise) + noise;
        [R_theta,R_rho] = cart2pol(real(Received_symbols),imag(Received_symbols)); 

        
        Last_symbols = zeros(1,N);
        for w=1:N
            if (R_theta(1,w) > pi)
                 Last_symbols(1,w) = R_rho(1,w)+1j*(R_theta(1,w)-2*pi);
            elseif (R_theta(1,w) < -pi)
                 Last_symbols(1,w) = R_rho(1,w)+1j*(R_theta(1,w)+2*pi);
            else
                 Last_symbols(1,w) = R_rho(1,w)+1j*(R_theta(1,w));
            end
        end
        
        
        %QAM(M(m)) part
        s_x_qam = real(transmitted_symbols_qam);
        s_y_qam = imag(transmitted_symbols_qam);

        [s_theta_qam,s_rho_qam] = cart2pol(s_x_qam,s_y_qam);
       

        const_symbols_qam = s_rho_qam+1j*s_theta_qam;
        
        received_symbols_qam = transmitted_symbols_qam.*exp(1j*phase_noise) + noise;
        [r_theta_qam,r_rho_qam] = cart2pol(real(received_symbols_qam),imag(received_symbols_qam));
        
        last_symbols_qam = zeros(1,N);
        for w=1:N
            if (r_theta_qam(1,w) > pi)
                 last_symbols_qam(1,w) = r_rho_qam(1,w)+1j*(r_theta_qam(1,w)-2*pi);
            elseif (r_theta_qam(1,w) < -pi)
                 last_symbols_qam(1,w) = r_rho_qam(1,w)+1j*(r_theta_qam(1,w)+2*pi);
            else
                 last_symbols_qam(1,w) = r_rho_qam(1,w)+1j*(r_theta_qam(1,w));
            end
        end
       
        
        %HQAM(M(m)) part
        s_x_hqam = real(transmitted_symbols_hqam);
        s_y_hqam = imag(transmitted_symbols_hqam);

        [s_theta_hqam,s_rho_hqam] = cart2pol(s_x_hqam,s_y_hqam);

        const_symbols_hqam = s_rho_hqam+1j*s_theta_hqam;
        
        received_symbols_hqam = transmitted_symbols_hqam.*exp(1j*phase_noise) + noise;
        [r_theta_hqam,r_rho_hqam] = cart2pol(real(received_symbols_hqam),imag(received_symbols_hqam));


        last_symbols_hqam = zeros(1,N);
        for w=1:N
            if (r_theta_hqam(1,w) > pi)
                 last_symbols_hqam(1,w) = r_rho_hqam(1,w)+1j*(r_theta_hqam(1,w)-2*pi);
            elseif (r_theta_hqam(1,w) < -pi)
                 last_symbols_hqam(1,w) = r_rho_hqam(1,w)+1j*(r_theta_hqam(1,w)+2*pi);
            else
                 last_symbols_hqam(1,w) = r_rho_hqam(1,w)+1j*(r_theta_hqam(1,w));
            end
        end
        
        % Set axis limits 
        
      
        x_axis_limit_1 = [-sqrt(2)*(x_max+sigma), sqrt(2)*(x_max+sigma)];
        y_axis_limit_1 = [-sqrt(2)*(x_max+sigma), sqrt(2)*(x_max+sigma)];

        image = plot(real(received_symbols) , imag(received_symbols) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_1);
        ylim(y_axis_limit_1);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images{((i-1)*num_images/length(No))+(m-1)*num_images+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        

         % Set axis limits 
        
        
        x_axis_limit_2 = [0, sqrt(2)*(x_max+sigma)];
        y_axis_limit_2 = [-pi , pi];

        image_v2 = plot(real(last_symbols) , imag(last_symbols) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_2);
        ylim(y_axis_limit_2);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        images_v2{((i-1)*num_images/length(No))+(m-1)*num_images+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        
        
        % Set axis limits 
       
        x_axis_limit_3 = [-sqrt(2)*(x_max+sigma), sqrt(2)*(x_max+sigma)];
        y_axis_limit_3 = [-sqrt(2)*(x_max+sigma), sqrt(2)*(x_max+sigma)];

        Image = plot(real(Received_symbols) , imag(Received_symbols) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_3);
        ylim(y_axis_limit_3);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        % Set axis limits using xlim and ylim
        images_q{((i-1)*num_images/length(No))+(m-1)*num_images+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        

         % Set axis limits 
        
        x_axis_limit_4 = [0, sqrt(2)*(x_max+sigma)];
        y_axis_limit_4 = [-pi , pi];

        Image_v2 = plot(real(Last_symbols) , imag(Last_symbols) , '.');
        axis off; % Turn off axis labels
        xlim(x_axis_limit_4);
        ylim(y_axis_limit_4);
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
        % Set axis limits using xlim and ylim
        images_q_v2{((i-1)*num_images/length(No))+(m-1)*num_images+k} = im2double(im2gray_new(frame2im(getframe)));
         
        % Set axis limits 
        
        
        x_axis_limit_7 = [-sqrt(2)*(x_max+sigma), sqrt(2)*(x_max+sigma)];
        y_axis_limit_7 = [-sqrt(2)*(x_max+sigma), sqrt(2)*(x_max+sigma)];
        image_qam = plot(real(received_symbols_qam) , imag(received_symbols_qam) , '.');
        axis off; % Turn off axis labels
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
         % Set axis limits using xlim and ylim
        xlim(x_axis_limit_7);
        ylim(y_axis_limit_7);
        images_qam{((i-1)*num_images/length(No))+(m-1)*num_images+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image

         % Set axis limits 
        
        x_axis_limit_8 = [0, sqrt(2)*(x_max+sigma)];
        y_axis_limit_8 = [-pi , pi];
        image_qam_v2 = plot(real(last_symbols_qam) , imag(last_symbols_qam) , '.');
        axis off; % Turn off axis labels
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
         % Set axis limits using xlim and ylim
        xlim(x_axis_limit_8);
        ylim(y_axis_limit_8);
        images_qam_v2{((i-1)*num_images/length(No))+(m-1)*num_images+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        
         % Set axis limits 
        
        
        x_axis_limit_11 = [-sqrt(2)*(x_max+sigma), sqrt(2)*(x_max+sigma)];
        y_axis_limit_11 = [-sqrt(2)*(x_max+sigma), sqrt(2)*(x_max+sigma)];
        image_hqam = plot(real(received_symbols_hqam) , imag(received_symbols_hqam) , '.');
        axis off; % Turn off axis labels
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
         % Set axis limits using xlim and ylim
        xlim(x_axis_limit_11);
        ylim(y_axis_limit_11);
        images_hqam{((i-1)*num_images/length(No))+(m-1)*num_images+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image

          % Set axis limits 
        
        x_axis_limit_12 = [0, sqrt(2)*(x_max+sigma)];
        y_axis_limit_12 = [-pi , pi];
        image_hqam_v2 = plot(real(last_symbols_hqam) , imag(last_symbols_hqam) , '.');
        axis off; % Turn off axis labels
        set(gca,'position',[0 0 1 1],'units','normalized') % Set plot position
        set(gcf,'position',[100 100 256 256]) % Set plot size
         % Set axis limits using xlim and ylim
        xlim(x_axis_limit_12);
        ylim(y_axis_limit_12);
        images_hqam_v2{((i-1)*num_images/length(No))+(m-1)*num_images+k} = im2double(im2gray_new(frame2im(getframe))); % Convert plot to image
        
    end
  end
end


%Create label matrix
labels_before = [];
for m = 1:length_M*4
    labels_before = [labels_before; m*ones(num_images, 1)];
end


num_classes = length_M*4;
labels_onehot_before = zeros(length(labels_before), num_classes);

for i = 1:length(labels_before)
    labels_onehot_before(i, labels_before(i)) = 1;
end


labels_onehot_total_before = labels_onehot_before;

%Create label matrix
labels_after = [];
for m = 1:length_M*4
    labels_after = [labels_after; m*ones(num_images, 1)];
end


num_classes = length_M*4;
labels_onehot_after = zeros(length(labels_after), num_classes);

for i = 1:length(labels_after)
    labels_onehot_after(i, labels_after(i)) = 1;
end

labels_onehot_total_after = labels_onehot_after;


%Convert the cell array to a regular array
Images = cell2mat(images);
Images_v2 = cell2mat(images_v2);

Images_q = cell2mat(images_q);
Images_q_v2 = cell2mat(images_q_v2);


Images_qam = cell2mat(images_qam);
Images_qam_v2 = cell2mat(images_qam_v2);

Images_hqam = cell2mat(images_hqam);
Images_hqam_v2 = cell2mat(images_hqam_v2);

Images_total_before = [Images; Images_q; Images_qam; Images_hqam;];
Images_total_after = [Images_v2; Images_q_v2; Images_qam_v2; Images_hqam_v2;];

image_order = zeros(1, num_images*length_M);
idx = 1;
for m=1:length_M
    for snr=min(SNR_dB):step_snr:max(SNR_dB)
        for i = 1:(num_images/length(No))
            image_order(idx) = snr;
            idx = idx + 1;
        end
    end
end
Image_order = transpose(image_order);

image_order_v2 = zeros(1, num_images*length_M);
idx = 1;
for m=1:length_M
    for snr=min(SNR_dB):step_snr:max(SNR_dB)
        for i = 1:(num_images/length(No))
            image_order_v2(idx) = snr;
            idx = idx + 1;
        end
    end
end
Image_order_v2 = transpose(image_order_v2);

image_order_q = zeros(1, num_images*length_M);
idx = 1;
for m=1:length_M
    for snr=min(SNR_dB):step_snr:max(SNR_dB)
        for i = 1:(num_images/length(No))
            image_order_q(idx) = snr;
            idx = idx + 1;
        end
    end
end
Image_order_q = transpose(image_order_q);

image_order_q_v2 = zeros(1, num_images*length_M);
idx = 1;
for m=1:length_M
    for snr=min(SNR_dB):step_snr:max(SNR_dB)
        for i = 1:(num_images/length(No))
            image_order_q_v2(idx) = snr;
            idx = idx + 1;
        end
    end
end
Image_order_q_v2 = transpose(image_order_q_v2);


image_order_qam = zeros(1, num_images*length_M);
idx = 1;
for m=1:length_M
    for snr=min(SNR_dB):step_snr:max(SNR_dB)
        for i = 1:(num_images/length(No))
            image_order_qam(idx) = snr;
            idx = idx + 1;
        end
    end
end
Image_order_qam = transpose(image_order_qam);

image_order_qam_v2 = zeros(1, num_images*length_M);
idx = 1;
for m=1:length_M
    for snr=min(SNR_dB):step_snr:max(SNR_dB)
        for i = 1:(num_images/length(No))
            image_order_qam_v2(idx) = snr;
            idx = idx + 1;
        end
    end
end
Image_order_qam_v2 = transpose(image_order_qam_v2);


image_order_hqam = zeros(1, num_images*length_M);
idx = 1;
for m=1:length_M
    for snr=min(SNR_dB):step_snr:max(SNR_dB)
        for i = 1:(num_images/length(No))
            image_order_hqam(idx) = snr;
            idx = idx + 1;
        end
    end
end
Image_order_hqam = transpose(image_order_hqam);

image_order_hqam_v2 = zeros(1, num_images*length_M);
idx = 1;
for m=1:length_M
    for snr=min(SNR_dB):step_snr:max(SNR_dB)
        for i = 1:(num_images/length(No))
            image_order_hqam_v2(idx) = snr;
            idx = idx + 1;
        end
    end
end
Image_order_hqam_v2 = transpose(image_order_hqam_v2);


Image_order_total_before = [Image_order ; Image_order_q;  Image_order_qam ; Image_order_hqam ;];
Image_order_total_after = [Image_order_v2 ; Image_order_q_v2; Image_order_qam_v2;  Image_order_hqam_v2;];

save('C:\Users\admin\Documents\MATLAB\ProjectEvgen\AMC_codes\apsk_qam_hqam_test.mat' ,'Images_total_before' ,'Images_total_after','labels_onehot_total_before' , 'labels_onehot_total_after','Image_order_total_before','Image_order_total_after')