%----------------------------------------------------------------
% This code can be used to generate Figure 1_5
%----------------------------------------------------------------
% This file requires the following files to be present in the same
% directory:
%
% tot30.mat
clear all
close all

load tot30; % load the measured scattered field

% DEFINITION OF PARAMETERS
f=linspace(6,18,251)*1e9; %Form frequency vector 
BW=6e9; % Select the frequency window size
d=2e-9; %Select the time delay 

% DISPLAY THE FIELD ÝN JFT PLANE
[B,T,F]=stft(tot30,f,BW,50,d);
xlabel('--->Time (nsec)'); ylabel('--> Freq. (GHz)');
colorbar;
colormap(1-gray)
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
axis tight; 
xlabel('Time [ns]');
ylabel('Frequency [GHz]');

