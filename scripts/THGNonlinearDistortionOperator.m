%This script takes the "DASH" corrected transmission matrix and transforms
%the data into the nonlinear distortion operator (NDO) form. In this form, 
%the NDO aberration correction method is implemented from the first
%singular value of the SVD of the NDO, as described in the main text of the
%manuscript "Synthetic aperture holographic third harmonic generation
%microscopy".
%The script in its current form runs section-by-setion. -YRF March 2024

%% Shift DASH corrected Transmission matrix to a "descanned" space
%this is the x_o -> x_m domain shift described in main text

D_shift = zeros(501*501,size(D_DASH,2));
x_shift = zeros(size(D_DASH,2));
y_shift = zeros(size(D_DASH,2));
ind = 1;
%center position for each measurement varies slightly due to how data was
%originally cropped from the raw camera pixels (1024x1024) to the cropped
%data (501x501).
Center_shift_positionx = 158; %158 for bone, 181 for mouse tail, 106 for MoS2
Center_shift_positiony = 94; %94 for bone, 141 for mouse tail, 95 for Mos2
for i = 1:26;
    for ii = 1:26;
    D_temp = reshape(D_DASH(:,ind),[501,501]);
    D_temp = circshift(D_temp,10*(ii-1)-Center_shift_positionx,2); %pixel shift of 10 pixels corresponds to sample position shift of 2.76 microns
    D_temp = circshift(D_temp,10*(i-1)-Center_shift_positiony,1); 
    x_shift(ind) = 10*(ii-1); %to keep track of shift position
    y_shift(ind) = 10*(i-1);
    %figure(1);imagesc(abs(D_temp));colormap turbo;daspect([1,1,1])
    D_shift(:,ind) = D_temp(:);
    ind = ind+1;
    end
end

%% Fourier Transform of each hologram, x_m -> u_o , create NDO


Distortion_matrix = zeros(size(D_shift,1),size(D_shift,2));
for i = 1:size(D_shift,2);
    temp = reshape(D_shift(:,i),[501,501]);
    temp = fftshift(fft2(ifftshift(temp)));
    %figure(2);imagesc(abs(temp));colormap turbo;daspect([1,1,1])
    Distortion_matrix(:,i)=temp(:);
end
%% Take SVD of NDO


%[U,S,V] = svd(Distortion_matrix,'econ');
[U,~,~] = svd(Distortion_matrix,'econ'); %only need the left singular value

%% Image of first singular value 

figure(3);imagesc(reshape(angle(U(:,1)),[501,501]));colormap turbo;daspect([1,1,1])

%% Creating variable for first singular value 

FirstSingVal = reshape(U(:,1),[501,501]);
%% Correcting aberrations with first singular value
%taking the inverse of the first singular value and applying this to each
%column of the NDO, with a regularization parameter, .  
    
Distortion_matrix_corrected = zeros(501*501,size(Distortion_matrix,2));
regul_param = 0.2; %0.015 for mouse tail, 0.2 for developing bone, 0.15 for MoS2

for i = 1:size(Distortion_matrix,2);
    temp = reshape(Distortion_matrix(:,i),[501,501])./(FirstSingVal+regul_param);
    %figure(4);imagesc(abs(temp))
    Distortion_matrix_corrected(:,i) = temp(:);
end


%% Taking the NDO corrected data and transforming it back to spatial domain

Distortion_matrix_corrected_FT = zeros(501*501,size(Distortion_matrix,2));

for i = 1:size(Distortion_matrix,2);
    temp = reshape(Distortion_matrix_corrected(:,i),[501,501]);
    temp = fftshift(ifft2(ifftshift(temp)));
    %figure(5);imagesc(abs(temp))
    Distortion_matrix_corrected_FT(:,i)=temp(:);
end

%% Shifting the data back to the coordinates of the Transmission matrix

Distortion_Final = zeros(501*501,size(Distortion_matrix,2));
ind = 1;
for i = 1:26;
    for ii = 1:26;
    D_temp = reshape(Distortion_matrix_corrected_FT(:,ind),[501,501]);
    D_temp = circshift(D_temp,-10*(ii-1)+Center_shift_positionx,2);
    D_temp = circshift(D_temp,-10*(i-1)+Center_shift_positiony,1);
    %x_shift(ind) = 10*(ii-1);
    %y_shift(ind) = 10*(i-1);
    %figure(6);imagesc(abs(D_temp))
    Distortion_Final(:,ind) = D_temp(:);
    ind = ind+1;
    end
end
figure();imagesc(abs(reshape(sum(Distortion_Final,2),[501,501])));daspect([1,1,1])


%%
x = 1:501;
y = x';
[X,Y] = meshgrid(x,y);

%% Create spatial circular filter function to each hologram in the final corrected transmission matrix
% Note: here a simple circular function is created and centers the circular
% filter to the maximum intensity of each reconstructed hologram, an
% alternative method is to use the circshift function to move the circular
% filter to the hologram location.

    ind = 1;
    Distortion_Final_filtered = zeros(501*501,size(D_DASH,2));
    for j = 1:sqrt(size(D_DASH,2));
        for jj = 1:sqrt(size(D_DASH,2));
          
            maximum = max(max(abs(reshape(Distortion_Final(:,ind),[501,501]))));
            [yy,xx]=find(abs(reshape(Distortion_Final(:,ind),[501,501]))==maximum);          
            filter = exp(-(((X-xx)/(50/2)).^2 + ((Y-yy)/(50/2)).^2).^4);
            temp = reshape((Distortion_Final(:,ind)),[501,501]).*filter;
            Distortion_Final_filtered(:,ind) = temp(:); 
            ind = ind + 1;

        end
    end
figure();imagesc(abs(reshape(sum(Distortion_Final_filtered,2),[501,501])));
%% Truncated SVD for further denoising

[U,S,V] = svd(Distortion_Final_filtered,'econ');
%% Set truncation point

k = 250; % Number of singular values to retain set by calculating variance 
         % of final image as a function of k and selecting point where
         % variance changes minimally as a function of k. 
S_truncated = S;
S_truncated(k+1:end, k+1:end) = 0;

Distortion_Final_filtered_denoised = U * S_truncated * V';
figure();imagesc(abs(reshape(sum(Distortion_Final_filtered_denoised,2),[501,501])));daspect([1,1,1])