clearvars; close all; clc;

fs = 250;           % Sample frequency
N = 3000;           % Data length
t = (-N:N)/fs;      % Time scale
s = 0.01:0.1:2;     % Values of scaling factor

% Maxican hat wavelet
%...

wavelt = -0.8671 .* exp(-0.5*t.^2).*(t.^2-1);



figure;
plot(t,wavelt,'b');
xlabel("t (s)");
ylabel("\psi(t)");
title("Mother Wavelet");
grid on;

k = 1;
for i = 1:5
    figure;
    for j = 1:4
        subplot(1,4,j);
        s_wavelt = (1/sqrt(s(k)))* (-0.8671 .* exp(-0.5*(t/s(k)).^2).*((t/s(k)).^2-1));
        plot(t, s_wavelt);
        ylim([-4,8.5]);
        title_string = sprintf("Scaling factor = %.2f",s(k));
        title(title_string);
        grid on; 
        disp(sprintf("When scaling factor is %.2f, mean = %.2f , energy = %.2f",s(k),mean(s_wavelt),round(trapz(t, wavelt.^2),2)));
        k=k+1;
    end 
end







% Generating spectra of wavelets
Fwavelt = fft(wavelt)/length(wavelt);
hz = linspace(0,fs/2,floor(length(wavelt)/2)+1);
figure;
plot(hz,2*abs(Fwavelt(1:length(hz))));
xlabel('Frequency (Hz)'), ylabel('Amplitude');
title("Spectra of the Mother Wavelet")

k = 1;
for i = 1:5
    figure;
    for j = 1:4
        subplot(1,4,j);
        s_wavelt = (1/sqrt(s(k)))* (-0.8671 .* exp(-0.5*(t/s(k)).^2).*((t/s(k)).^2-1));

        Fwavelt = fft(s_wavelt)/length(s_wavelt);
        hz = linspace(0,fs/2,floor(length(s_wavelt)/2)+1);
        
        plot(hz,2*abs(Fwavelt(1:length(hz))));
        xlabel('Frequency (Hz)'), ylabel('Amplitude');
        ylim([0,0.2]);
        title_string = sprintf("Scaling factor = %.2f",s(k));
        title(title_string);
        grid on; 
        k=k+1;
    end 
end


% % Ploting the spectrogram
% h = pcolor();
% set(h, 'EdgeColor', 'none');
% colormap jet
% xlabel('Time (s)')
% ylabel('Scale')
% title('Spectrogram')