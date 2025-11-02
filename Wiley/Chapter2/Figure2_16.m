%----------------------------------------------------------------
% This code can be used to generate Figure 2.16
%----------------------------------------------------------------
clear all
close all

t=0:0.01e-9:50e-9; % choose time vector
N=length(t);
pulse(201:300)=ones(1,100); % form rectangular pulse
pulse(N)=0; 

%  Frequency domain equivalent
dt=t(2)-t(1); % time resolution
df=1/(N*dt); % frequency resolution
f=0:df:df*(N-1); % set frequency vector
pulseF=fft(pulse)/N; % frequency domain signal

%---Figure 2.16(a)------------------------------------------------
plot(t(1:501)*1e9,pulse(1:501),'k','LineWidth',2);
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
xlabel('Time [ns]'); ylabel('Amplitude [V]'); 
axis([0 5 0 1.1]); 

%---Figure 2.16(b)------------------------------------------------
plot(f(1:750)/1e6,20*log10(abs(pulseF(1:750))),'k','LineWidth',2);
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
axis([0 f(750)/1e6 -85 -20]);
xlabel('Frequency [MHz]'); ylabel('Amplitude [dB]'); 

