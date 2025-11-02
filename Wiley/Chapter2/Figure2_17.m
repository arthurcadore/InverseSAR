%----------------------------------------------------------------
% This code can be used to generate Figure 2.17
%----------------------------------------------------------------
clear all
close all

t=-25e-9:0.01e-9:25e-9; % choose time vector
N=length(t); 
sine(2451:2551)=-sin(2*pi*1e9*(t(2451:2551)));% form sine pulse
sine(N)=0; 

%  Frequency domain equivalent
dt=t(2)-t(1); % time resolution
df=1/(N*dt); % frequency resolution
f=0:df:df*(N-1); % set frequency vector
sineF=fft(sine)/N; % frequency domain signal

%---Figure 2.17(a)------------------------------------------------
plot(t(2251:2751)*1e9,sine(2251:2751),'k','LineWidth',2);
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
xlabel('Time [ns]'); ylabel('Amplitude [V]'); 
axis([-2.5 2.5 -1.1 1.1]); 

%---Figure 2.17(b)------------------------------------------------
plot(f(1:750)/1e6,20*log10(abs(sineF(1:750))),'k','LineWidth',2);
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
axis([0 f(750)/1e6 -85 -20]);
xlabel('Freq [ MHz]'); ylabel('Amplitude [dB]'); 

