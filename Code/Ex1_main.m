clc
clear all
close all

% Set params
SDS = 3; % [cm]
readtable("DPFperTissue.txt")
tissues = readtable("DPFperTissue.txt").Tissue;
tissueType = tissues{1};

% Run function on the first file
firstFile = '.\FN_032_V1_Postdose1_Nback.mat';
[ dHbR , dHbO, fig ] = CalcNIRS(firstFile, SDS, tissueType, [1, 2]);
fig.Name = tissueType;
fig.NumberTitle = 'off';       

% Run function on the second file
secondFile = '.\FN_031_V2_Postdose2_Nback.mat';
[ ~, ~, fig ] = CalcNIRS(secondFile, SDS, tissueType, [1, 2]);
fig.Name = tissueType;
fig.NumberTitle = 'off';       

% Calc the Fourier transform of the first channel in first file
FTdHbR = (fft(dHbR(:, 1)));
FTdHbO = (fft(dHbO(:, 1)));
% Calc frequency axis
load(firstFile);
N = length(t);
dt = mean(diff(t));
Fs = 1/dt;
freqs = (0:N-1)*(Fs/N);   % Frequency axis in Hz
% take the positive frequencies
FTdHbR = FTdHbR(1:N/2);
FTdHbO = FTdHbO(1:N/2);
freqs = freqs(1:N/2);

% Plot the Fourier transform magnitude
fig = figure;
hold;
plot(freqs, abs(FTdHbR), 'b'); 
plot(freqs, abs(FTdHbO), 'r'); 
title('Fourier transform of the first channel');
legend(["dHbR", "dHbO"]);
xlabel("Frequency (Hrz)");
ylabel("Fourier transform magnitude (t*mol/L)");

% Calc SNR based on the following rules:
% - Noise = average of the signal strength in Fourier domain at Frequencies above 2.5Hz
% - Signal = signal strength at heart beat frequency
% I assumed heart beat frequency must be between 40 and 150 (bpm) and the
% real heart beat is the peak in this frequncies range
hbfZone = find(all([40 <= freqs*60 ; freqs*60 <= 150])); 
noiseZone = find(freqs >= 2.5);
SNRdHbR = max(abs(FTdHbR(hbfZone))) / mean(abs(FTdHbR(noiseZone)));
SNRdHbO = max(abs(FTdHbO(hbfZone))) / mean(abs(FTdHbO(noiseZone)));
disp(['The SNR of dHbR is:  ', num2str(SNRdHbR)]);
disp(['The SNR of dHbO is:  ', num2str(SNRdHbO)]);


% Now check the effect of changing tissueType parameter
for i=1:4
    tissueType = tissues{i};
    % Run function 
    [ ~ , ~, fig ] = CalcNIRS(firstFile, SDS, tissueType, [1, 2]);
    fig.Name = tissueType;
end
