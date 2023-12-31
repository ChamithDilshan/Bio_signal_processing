function Hd = Kaiser_high_pass_application
%KAISER_HIGH_PASS_APPLICATION Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.14 and Signal Processing Toolbox 9.2.
% Generated on: 19-Sep-2023 14:10:39

% FIR Window Highpass filter designed using the FIR1 function.

% All frequency values are in Hz.
Fs = 500;  % Sampling Frequency

N    = 302;      % Order
Fc   = 5;        % Cutoff Frequency
flag = 'scale';  % Sampling Flag
Beta = 5.6533;   % Window Parameter

% Create the window vector for the design algorithm.
win = kaiser(N+1, Beta);

% Calculate the coefficients using the FIR1 function.
b  = fir1(N, Fc/(Fs/2), 'high', win, flag);
Hd = dfilt.dffir(b);

% [EOF]
