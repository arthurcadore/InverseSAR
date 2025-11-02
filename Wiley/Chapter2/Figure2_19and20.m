%----------------------------------------------------------------
% This code can be used to generate Figure 2.19 and Figure 2.20
%----------------------------------------------------------------
clear all
close all

fo=1e6;  % choose base frequency
t=0:1e-9:4e-6; tt=0:1e-9:17e-6;  % choose time vector
k=1.0e12;  % choose chirp rate
sinep=sin(2*pi*fo*t);sinep(12001:16001)=sinep; % form CW pulse
sinep(17001)=0;
m=sin(2*pi*(fo+k*t/2).*t); % form LFM pulse
m(12001:16001)=m; m(17001)=0;

%---Figure 2.19(a)------------------------------------------------
plot(tt*1e6,sinep,'k','LineWidth',2);
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
xlabel('Time [\mus]'); ylabel('Amplitude [V]'); 
axis([0 17 -1.1 1.1])

%---Figure 2.19(b)------------------------------------------------
plot(tt*1e6,m,'k','LineWidth',2);
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
xlabel('Time [\mus]'); ylabel('Amplitude [V]'); 
axis([0 17 -1.1 1.1])

%  Frequency domain equivalent
df=1/170e-6; % frequency resolution
f=0:df:df*170000; % set frequency vector
sinep=sin(2*pi*fo*t);sinep(170001)=0;
m=sin(2*pi*(fo+k*t/2).*t); m(170001)=0;

fsinep=fft(sinep)/length(t); % spectrum of CW pulse
fm=fft(m)/length(t);% spectrum of LFM pulse

%---Figure 2.20(a)------------------------------------------------
plot(f(1:2000)/1e6,20*log10(abs(fsinep(1:2000))),'k','LineWidth',2);
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
xlabel('Frequency [MHz]'); ylabel('Amplitude [dB]'); 
text(8,-10,'Single tone pulse')
axis([0 12 -40 -5])

%---Figure 2.20(b)------------------------------------------------
plot(f(1:2000)/1e6,20*log10(abs(fm(1:2000))),'k','LineWidth',2);
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
xlabel('Frequency [MHz]'); ylabel('Amplitude [dB]'); 
text(8,-10,'Chirp pulse')
axis([0 12 -40 -5])
