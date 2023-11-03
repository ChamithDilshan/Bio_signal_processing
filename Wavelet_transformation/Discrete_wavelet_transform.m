%% Discrete Wavelet Transform

clearvars; close all; clc;

n = 0:1023;
fs = 512;
%% 
t = n/fs;

x1 = [2*sin(20*pi*(0:511)/fs)+sin(80*pi*(0:511)/fs) 0.5*sin(40*pi*(512:1023)/fs)+sin(60*pi*(512:1023)/fs)]; %Construction of x1

x2 = zeros(1,1024);  %Construction of x2
for i=0:1023
    if (0<=i) && (i<64)
        x2(i+1) = 1;
    elseif (192<=i) && (i<256)
        x2(i+1) = 2;
    elseif (256<=i) && (i<512)
        x2(i+1) = -1;
    elseif (512<=i) && (i<704)
        x2(i+1) = 3;
    elseif (704<=i) && (i<960)
        x2(i+1) = 1;
    else
        x2(i+1) = 0;
    end
end

rng(10); %Seed
y1 = awgn(x1,10,"measured");  %Adding noise to signals
y2 = awgn(x2,10,"measured");

%Plotting the signals
figure;
plot(t,x1);
hold on;
plot(t,y1);
legend("x_1[n]","y_1[n]");
title("x_1[n] and y_1[n]");

figure;
plot(t,x2);
hold on;
plot(t,y2);
legend("x_2[n]","y_2[n]");
title("x_2[n] and y_2[n]");

%% Morphology of Wavelets

figure;

subplot(1,2,1);
[phi, psi, x] = wavefun('haar', 10); 
plot(x, phi, 'b', x, psi, 'r');
legend('Scaling Function (phi)', 'Wavelet Function (psi)');
title('Haar Wavelet and Scaling Functions');

subplot(1,2,2)
[phi, psi, x] = wavefun('db9', 10);  
plot(x, phi, 'b', x, psi, 'r');
legend('Scaling Function (phi)', 'Wavelet Function (psi)');
title('Daubechies D9 Wavelet and Scaling Functions');

waveletAnalyzer

%% Decomposition using code  

% Decomposition of signal y1
[c_db9_y1, l_db9_y1] = wavedec(y1, 10, 'db9');  
[c_haar_y1, l_haar_y1] = wavedec(y1, 10, 'haar'); 

%Decomposition of signal y2
[c_db9_y2, l_db9_y2] = wavedec(y2, 10, 'db9');  
[c_haar_y2, l_haar_y2] = wavedec(y2, 10, 'haar');

%% Using inverse DWT to reconstruct the signal

reconstructed_db9_y1 = waverec(c_db9_y1, l_db9_y1, 'db9');
energy_db9_y1 = sum((y1 - reconstructed_db9_y1).^2);

reconstructed_haar_y1 = waverec(c_haar_y1, l_haar_y1, 'haar');
energy_haar_y1 = sum((y1 - reconstructed_haar_y1).^2);

reconstructed_db9_y2 = waverec(c_db9_y2, l_db9_y2, 'db9');
energy_db9_y2 = sum((y2 - reconstructed_db9_y2).^2);

reconstructed_haar_y2 = waverec(c_haar_y2, l_haar_y2, 'haar');
energy_haar_y2 = sum((y2 - reconstructed_haar_y2).^2);

figure;
subplot(2,2,1)
plot(t,y1);
hold on;
plot(t,reconstructed_db9_y1);
legend("y_1[n]","Reconstructed y_1 with db9",'Location', 'SouthEast');

subplot(2,2,2)
plot(t,y1);
hold on;
plot(t,reconstructed_haar_y1);
legend("y_1[n]","Reconstructed y_1 with haar",'Location', 'SouthEast');

subplot(2,2,3)
plot(t,y2);
hold on;
plot(t,reconstructed_db9_y2);
legend("y_2[n]","Reconstructed y_2 with db9",'Location', 'SouthEast');

subplot(2,2,4)
plot(t,y2);
hold on;
plot(t,reconstructed_haar_y2);
legend("y_2[n]","Reconstructed y_2 with haar",'Location', 'SouthEast');

%% Signal Denoising with DWT

clc;

% ------------------------------------------------For the signal y1 with db9------------------------------------------------

magnitude_db9_y1 = abs(c_db9_y1); %get magnitudes
[magnitude_db9_y1_sorted, idx] = sort(magnitude_db9_y1, 'descend'); % Sort the coefficients and their magnitudes in descending order
coefficients_db9_y1_sorted = c_db9_y1(idx);

threshold_value = 1;
thresholded_coefficients = c_db9_y1;
thresholded_coefficients(abs(c_db9_y1) < threshold_value) = 0;
reconstructed_signal_db9_y1 = waverec(thresholded_coefficients, l_db9_y1, 'db9');

fprintf("Signal y_1[n], db9 wavelet ----> RMSE = %.3f\n",mean((x1-reconstructed_signal_db9_y1).^2));

% ------------------------------------------------For the signal y2 with db9------------------------------------------------

magnitude_db9_y2 = abs(c_db9_y2); %get magnitudes
[magnitude_db9_y2_sorted, idx] = sort(magnitude_db9_y2, 'descend'); % Sort the coefficients and their magnitudes in descending order
coefficients_db9_y2_sorted = c_db9_y2(idx);

threshold_value = 1.7;
thresholded_coefficients = c_db9_y2;
thresholded_coefficients(abs(c_db9_y2) < threshold_value) = 0;
reconstructed_signal_db9_y2 = waverec(thresholded_coefficients, l_db9_y2, 'db9');

fprintf("Signal y_2[n], db9 wavelet ----> RMSE = %.3f\n",mean((x2-reconstructed_signal_db9_y2).^2));

% ------------------------------------------------For the signal y1 with haar------------------------------------------------

magnitude_haar_y1 = abs(c_haar_y1); %get magnitudes
[magnitude_haar_y1_sorted, idx] = sort(magnitude_haar_y1, 'descend'); % Sort the coefficients and their magnitudes in descending order
coefficients_haar_y1_sorted = c_haar_y1(idx);

threshold_value = 0.1;
thresholded_coefficients = c_haar_y1;
thresholded_coefficients(abs(c_haar_y1) < threshold_value) = 0;
reconstructed_signal_haar_y1 = waverec(thresholded_coefficients, l_haar_y1, 'haar');

fprintf("Signal y_1[n], haar wavelet ----> RMSE = %.3f\n",mean((x1-reconstructed_signal_haar_y1).^2));

% ------------------------------------------------For the signal y2 with haar------------------------------------------------

magnitude_haar_y2 = abs(c_haar_y2); %get magnitudes
[magnitude_haar_y2_sorted, idx] = sort(magnitude_haar_y2, 'descend'); % Sort the coefficients and their magnitudes in descending order
coefficients_haar_y2_sorted = c_haar_y2(idx);

threshold_value = 2.1;
thresholded_coefficients = c_haar_y2;
thresholded_coefficients(abs(c_haar_y2) < threshold_value) = 0;
reconstructed_signal_haar_y2 = waverec(thresholded_coefficients, l_haar_y2, 'haar');

fprintf("Signal y_2[n], haar wavelet ----> RMSE = %.3f\n",mean((x2-reconstructed_signal_haar_y2).^2));



%-----------Plotting stems for coefficients-----------
figure;

subplot(2, 2, 1);
stem(coefficients_db9_y1_sorted, 'Marker', 'none');
title('Magnitude of db9 Wavelet Coefficients for y1 (Descending Order)');
xlabel('Coefficient Index');
ylabel('Magnitude');

subplot(2, 2, 2);
stem(coefficients_db9_y2_sorted, 'Marker', 'none');
title('Magnitude of db9 Wavelet Coefficients for y2 (Descending Order)');
xlabel('Coefficient Index');
ylabel('Magnitude');

subplot(2, 2, 3);
stem(coefficients_haar_y1_sorted, 'Marker', 'none');
title('Magnitude of haar Wavelet Coefficients for y1 (Descending Order)');
xlabel('Coefficient Index');
ylabel('Magnitude');

subplot(2, 2, 4);
stem(coefficients_haar_y2_sorted, 'Marker', 'none');
title('Magnitude of haar Wavelet Coefficients for y2 (Descending Order)');
xlabel('Coefficient Index');
ylabel('Magnitude');

%-----------Plotting denoised signals-----------
figure;

subplot(2, 2, 1);
plot(x1);
hold on;
plot(reconstructed_signal_db9_y1);
title("y_1[n] denoising with db9");
legend("x_1[n]", "Denoised");

subplot(2, 2, 2);
plot(x2);
hold on;
plot(reconstructed_signal_db9_y2);
title("y_2[n] denoising with db9");
legend("x_2[n]", "Denoised");

subplot(2, 2, 3);
plot(x1);
hold on;
plot(reconstructed_signal_haar_y1);
title("y_1[n] denoising with haar");
legend("x_1[n]", "Denoised");

subplot(2, 2, 4);
plot(x2);
hold on;
plot(reconstructed_signal_haar_y2);
title("y_2[n] denoising with haar");
legend("x_2[n]", "Denoised");


%% Signal Compression with DWT

clearvars;

load("ECGsig.mat");

fs = 257;
t = (0:length(aVR)-1)/fs;

figure;
plot(t,aVR);
xlabel("Time");
ylabel("Amplitude");
title("Given _aV_R ECG");

[c_db9, l_db9] = wavedec(aVR,5, 'db9');  %decomposition with 5 layers
[c_haar, l_haar] = wavedec(aVR,5, 'haar');

magnitude_db9 = abs(c_db9); %get magnitudes
[magnitude_db9_sorted, idx_db9] = sort(magnitude_db9, 'descend'); % Sort the coefficients and their magnitudes in descending order
coefficients_db9_sorted = c_db9(idx_db9);

magnitude_haar = abs(c_haar); 
[magnitude_haar_sorted, idx_haar] = sort(magnitude_haar, 'descend'); 
coefficients_haar_sorted = c_haar(idx_haar);

%Calculate the energy for db9
energy_db9 = zeros(size(coefficients_db9_sorted));
for coefs = 1:length(coefficients_db9_sorted)
    if coefs == 1
        energy_db9(coefs) = coefficients_db9_sorted(coefs)^2;
    else
        energy_db9(coefs) = coefficients_db9_sorted(coefs)^2 + energy_db9(coefs-1);
    end
end
energy_db9 = energy_db9 / max(energy_db9);

index_99_db9 = find(energy_db9 >= 0.99,1); %find the number of coefficients for 99% energy


%Calculate the energy for haar
energy_haar = zeros(size(coefficients_haar_sorted));
for coefs = 1:length(coefficients_haar_sorted)
    if coefs == 1
        energy_haar(coefs) = coefficients_haar_sorted(coefs)^2;
    else
        energy_haar(coefs) = coefficients_haar_sorted(coefs)^2 + energy_haar(coefs-1);
    end
end
energy_haar = energy_haar / max(energy_haar);

index_99_haar = find(energy_haar >= 0.99,1); %find the number of coefficients for 99% energy

figure;
subplot(1,2,1);
plot(energy_db9); hold on;
plot(index_99_db9,energy_db9(index_99_db9),'xr');
ylim([0 1.1]); xlabel("Number of coefficients"); ylabel("Energy percentage"); title("Energy vs number of coefficients");
text(index_99_db9, energy_db9(index_99_db9), sprintf('(%d,%.2f)', index_99_db9,energy_db9(index_99_db9)));

subplot(1,2,2);
plot(energy_haar); hold on;
plot(index_99_haar,energy_haar(index_99_haar),'xr');
ylim([0 1.1]); xlabel("Number of coefficients"); ylabel("Energy percentage"); title("Energy vs number of coefficients");
text(index_99_haar, energy_haar(index_99_haar), sprintf('(%d,%.2f)', index_99_haar,energy_haar(index_99_haar)));

% Reconstruction after compression
new_coefficients_db9 = zeros(size(coefficients_db9_sorted));
for i = 1:index_99_db9
    new_coefficients_db9(idx_db9(i)) = coefficients_db9_sorted(i);
end
compressed_db9 = waverec(new_coefficients_db9, l_db9, 'db9');

% Reconstruction after compression
new_coefficients_haar = zeros(size(coefficients_haar_sorted));
for i = 1:index_99_haar
    new_coefficients_haar(idx_haar(i)) = coefficients_haar_sorted(i);
end
compressed_haar = waverec(new_coefficients_haar, l_haar, 'haar');

figure;
plot(aVR);
hold on;
plot(compressed_haar);
hold on;
plot(compressed_db9);
legend("_aV_R","Compressed using haar","Compressed using db9");

fprintf("With haar wavelet ----> RMSE = %.7f , compression ratio = %d:1\n",mean((aVR-compressed_haar).^2),round(length(aVR)/index_99_haar));
fprintf("With db9 wavelet ----> RMSE = %.7f , compression ratio = %d:1\n",mean((aVR-compressed_db9).^2),round(length(aVR)/index_99_db9));
