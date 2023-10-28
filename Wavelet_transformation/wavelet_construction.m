clearvars; close all; clc;

fs = 250;           % Sample frequency
N = 3000;           % Data length
t = (-N:N)/fs;      % Time scale
s = 0.01:0.1:2;     % Values of scaling factor

% Maxican hat wavelet
%...

wavelt = 0.8671 .* exp(-0.5*t.^2).*(t.^2-1);



figure;
plot(t,wavelt,'b');
xlabel("t (s)");
ylabel("\psi(t)");

hold on;

legend_text = cell(length(s)+1, 1);
legend_text{1} = 'Original';

for i = 1:length(s)
    s_wavelt = (1/sqrt(s(i)))*wavelt;
    plot(t, s_wavelt);
    legend_text{i+1} = ['Scale Factor ' num2str(s(i))];
    hold on;
end

legend(legend_text);




% % Generating spectra of wavelets
% Fwavelt = fft(wavelt)/length(wavelt);
% hz = linspace(0,fs/2,floor(length(wavelt)/2)+1);
% plot(hz,2*abs(Fwavelt(1:length(hz))))
% xlabel('Frequency (Hz)'), ylabel('Amplitude')

% % Ploting the spectrogram
% h = pcolor();
% set(h, 'EdgeColor', 'none');
% colormap jet
% xlabel('Time (s)')
% ylabel('Scale')
% title('Spectrogram')