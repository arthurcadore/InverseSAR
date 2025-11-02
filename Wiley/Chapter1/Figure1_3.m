%----------------------------------------------------------------
% This code can be used to generate Figure 1_3
%----------------------------------------------------------------
% This file requires the following files to be present in the same
% directory:
%
%prince.wav

clear all
close all

% Read the sound signal "prince.wav"
[y,Fs,bits] = wavread('prince.wav');
sound(y,Fs);
N = length(y);

t = 0:.8/(N-1):.8; %form time vector
df = 1/max(t);
f = 0:df:df*(length(t)-1);

% TIME FREQUENCY PLANE SIGNAL              
A=spectrogram(y,256,250,400,1e4); % Calculate the spectrogram
matplot(t,f, (abs(A)),30); % Display the signal in T-F domain
colormap(1-gray);         % Change the colormap to grayscale
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
xlabel('Time [s]');
ylabel('Frequency [Hz]');
title('signal in time-frequency plane');
