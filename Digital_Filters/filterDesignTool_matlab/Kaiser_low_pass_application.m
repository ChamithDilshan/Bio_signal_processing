function Hd = Kaiser_low_pass_application

Fs = 500;  % Sampling Frequency

N    = 61;       % Order
Fc   = 125;      % Cutoff Frequency
flag = 'scale';  % Sampling Flag
Beta = 5.6533;   % Window Parameter

% Create the window vector for the design algorithm.
win = kaiser(N+1, Beta);

% Calculate the coefficients using the FIR1 function.
b  = fir1(N, Fc/(Fs/2), 'low', win, flag);


freqz(b,1,1024,500);

Hd = dfilt.dffir(b);
% [EOF]


