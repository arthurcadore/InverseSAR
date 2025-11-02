%----------------------------------------------------------------
% This code can be used to generate Figure 1_2
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

t=0:.8/(N-1):.8; %form time vector

% FREQUENCY DOMAIN SIGNAL
Y=fft(y)/N;
% Calculate the spectrum of the signal
df=1/(max(t)-min(t));           % Find the resolution in frequency   
f=0:df:df*(length(t)-1);        % Form the frequency vector
plot(f(1:2:N),abs(Y(1:(N+1)/2)),'k')  %downsample for plotting
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
axis tight; 
xlabel('Frequency [Hz]');
ylabel('Amplitude');
title('frequency domain signal');
