%% Discrete Wavelet Transform

clearvars; close all; clc;

n = 0:1023;
fs = 512;
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

%% 
