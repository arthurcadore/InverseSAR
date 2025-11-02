%----------------------------------------------------------------
% This code can be used to generate Figure 2.18
%----------------------------------------------------------------
clear all
close all

sigma=1e-10; % set sigma
t=-25e-9:0.01e-9:25e-9;  % choose time vector
N=length(t); 
mex(2451:2551)=1/sqrt(2*pi)/sigma^3*(1-t(2451:2551).^2/sigma^2)...
    .*(exp(-t(2451:2551).^2/2/sigma^2)); % form wavelet
mex=mex/max(mex);mex(N)=0; 

%  Frequency domain equivalent
dt=t(2)-t(1); % time resolution
df=1/(N*dt); % frequency resolution
f=0:df:df*(N-1); % set frequency vector
mexF=fft(mex)/N; % frequency domain signal

%---Figure 2.18(a)------------------------------------------------
plot(t(2251:2751)*1e9,mex(2251:2751),'k','LineWidth',2);
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
xlabel('Time [ns]'); ylabel('Amplitude [V]'); 
axis([-2.5 2.5 -.5 1.1]); 

%---Figure 2.18(a)------------------------------------------------
plot(f(1:750)/1e6,20*log10(abs(mexF(1:750))),'k','LineWidth',2);
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
axis([0 f(750)/1e6 -180 -40]);
xlabel('Freq [ MHz]'); ylabel('Amplitude [dB]'); 
