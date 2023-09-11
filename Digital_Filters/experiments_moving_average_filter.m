close all;
load("ECG_template.mat");

fs = 500;  % Sampling frequency = 500 Hz
T = 1/fs;  % T= Sampling period

num_points = size(ECG_template, 2);
time_axis = linspace(0, T * (num_points-1), num_points);  %time axis from the frequency

figure;
plot(time_axis, ECG_template);
xlabel('Time (seconds)');
ylabel('Amplitude (mv)');
title('ECG Template');

%---------------Noisy ECG---------------------------------

snr = 5; %signal-noise ratio in decible
nECG = awgn(ECG_template,snr,'measured');

figure;
plot(time_axis, nECG);
xlabel('Time (seconds)');
ylabel('Amplitude (mv)');
title('ECG with noise');

%---------------PSD of Noisy ECG------------------------------------
nfft = num_points;
[pxx_nECG,f_nECG] = periodogram(nECG,window,nfft,fs);

figure;
plot(f_nECG, 10*log10(pxx_nECG)); % Plot in dB scale
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('Power Spectral Density (PSD) of Noisy ECG');

%----------Moving Average Filter--------------------------------

%------MA(3) filter implementation with a customised script------

ma3ECG_1 = moving_average_filter(3,nECG);

group_delay = 3/2;
compensated_time_axis = time_axis - group_delay*T;

figure;
plot(time_axis, nECG);
hold on;  
plot(time_axis, ma3ECG_1, 'g--');
hold on;  
plot(compensated_time_axis, ma3ECG_1, 'r');
xlabel('Time (seconds)');
ylabel('Amplitude (mv)');
title('nECG and scripted MA(3)');
legend('nECG', 'ma3ECG_1','compensated');
hold off;


%------Compare PSDs------

[pxx_ma3,f_ma3] = periodogram(ma3ECG_1,window,nfft,fs);

figure;
plot(f_nECG, 10*log10(pxx_nECG));
hold on;
plot(f_ma3, 10*log10(pxx_ma3),'r');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('Power Spectral Densities (PSDs)');
legend('nECG', 'ma3ECG_1');
hold off;


%----MA(3) filter implementation with the MATLAB built-in function-----
windowSize = 3; 
numerator = (1/windowSize)*ones(1,windowSize);
denominator = 1;

ma3ECG_2 = filter(numerator,denominator,nECG);

figure;
plot(time_axis, nECG);
hold on;  
plot(time_axis, ma3ECG_2, 'r');
xlabel('Time (seconds)');
ylabel('Amplitude (mv)');
title('nECG and builtin MA(3)');
legend('nECG', 'ma3ECG_2');
hold off;




