%----------------------------------------------------------------
% This code can be used to generate Figure 2.11
%----------------------------------------------------------------
clear all
close all

fo=100; % set the base  frequency
t=0:1e-7:4e-3; % choose time vector
k=3e6; % select chirp rate
m=sin(2*pi*(fo+k*t/2).*t); % time domain FMCW signal

plot(t*1e3,m,'k','LineWidth',2);
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
xlabel('Time [ms]'); ylabel('Amplitude [V]'); 
axis([0 4 -1.1 1.1])


