%----------------------------------------------------------------
% This code can be used to generate Figure 1_1
%----------------------------------------------------------------
% This file requires the following files to be present in the same
% directory:
%
%prince.wav

clear all
close all

% Read the sound signal "prince.wav"
[y,Fs,bits] = wavread('prince.wav');
sound(y,Fs); %play the sound
N=length(y);

% TIME DOMAIN SIGNAL
t=0:.8/(N-1):.8; %form time vector
plot(t,y,'k');  %downsample for plotting
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
axis tight; 
xlabel('Time [s]');
ylabel('Amplitude');
title('time domain signal');
