%----------------------------------------------------------------
% This code can be used to generate Figure 2.15
%----------------------------------------------------------------
clear all
close all

c=.3; % speed of light []
f=2:.002:22; % choose frequency vector
Ro=50; % choose range of target [m]
k=2*pi*f/c;

Es=1*exp(-j*2*k*Ro); % collected SFCW electric field

df=f(2)-f(1); % frequency resolution
N=length(f); % total stepped frequency points
dr=c/(2*N*df); % range resolution
R=0:dr:dr*(length(f)-1); %set the range vector

plot(R, 20*log10(abs(ifft(Es))),'k','LineWidth',2)
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
xlabel('Range [m]'); ylabel('Amplitude [dBsm]'); 
axis([0 max(R) -90 0])
