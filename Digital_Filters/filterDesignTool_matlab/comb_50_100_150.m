function Hd = comb_50_100_150

Fs = 500;  % Sampling Frequency

% Notch Filter 1: 50 Hz
Fnotch1 = 50;   % Notch Frequency 1
BW1 = 3;        % Bandwidth 1
Apass1 = 1;     % Bandwidth Attenuation 1
[b1, a1] = iirnotch(Fnotch1 / (Fs/2), BW1 / (Fs/2), Apass1);

% Notch Filter 2: 100 Hz
Fnotch2 = 100;  % Notch Frequency 2
BW2 = 3;        % Bandwidth 2
Apass2 = 1;     % Bandwidth Attenuation 2
[b2, a2] = iirnotch(Fnotch2 / (Fs/2), BW2 / (Fs/2), Apass2);

% Notch Filter 3: 150 Hz
Fnotch3 = 150;  % Notch Frequency 3
BW3 = 3;        % Bandwidth 3
Apass3 = 1;     % Bandwidth Attenuation 3
[b3, a3] = iirnotch(Fnotch3 / (Fs/2), BW3 / (Fs/2), Apass3);

% Combine the filters
b_combined = conv(conv(b1, b2), b3);
a_combined = conv(conv(a1, a2), a3);

freqz(b_combined, a_combined, 1024, Fs);
Hd = dfilt.df2(b_combined, a_combined);

% [EOF]
