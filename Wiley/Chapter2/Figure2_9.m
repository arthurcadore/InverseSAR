%----------------------------------------------------------------
% This code can be used to generate Figure 2.9
%----------------------------------------------------------------
%---Figure 2.9(a)------------------------------------------------
clear all
close all

fo=1e3; % set the frequency
t=-4e-3:1e-7:4e-3; % choose time vector
s=cos(2*pi*fo*t); %  time domain CW signal
plot(t*1e3,s,'k','LineWidth',2);
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
xlabel('Time [ms]');
ylabel('Amplitude [V]');
axis([-4 4 -1.2 1.2])

%---Figure 2.9(b)------------------------------------------------
N=length(t);
df=1/(t(N)-t(1)); % Find frequency resolution
f=-df*(N-1)/2:df:df*(N-1)/2; % set frequency vector

S=fft(s)/N; % frequency domain CW signal
plot(f,fftshift(abs(S)),'k','LineWidth',2);
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
xlabel('Frequency [Hz]');
ylabel('Amplitude [V]');
axis([-.8e4 .8e4 0 .6])


