%----------------------------------------------------------------
% This code can be used to generate Figure 1_11
%----------------------------------------------------------------
clear all
close all

% TIME DOMAIN SIGNAL
a=0:.1:1;
t=(0:10)*1e-3;
stem(t*1e3,a,'k','Linewidth',2);% Figure 1-11 (a)
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
xlabel('time [ms]'); ylabel('s[n]');axis([-0.2 10.2 0 1.2]);


% FREQUENCY DOMAIN SIGNAL
b=fft(a);
df=1./(t(11)-t(1)); f=(0:10)*df; ff=(-5:5)*df; 
figure;
stem(f,abs(b),'k','Linewidth',2); % Figure 1-11 (b)
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
xlabel('frequency [Hz]'); ylabel('S[k]');axis([-20 1020 0 6.5]);
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');

figure;
stem(ff,fftshift(abs(b)),'k','Linewidth',2);% Figure 1-11 (c)
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
xlabel('frequency [Hz]'); ylabel('S[k]');axis([-520 520 0 6.5]);
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');


