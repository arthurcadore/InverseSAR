%----------------------------------------------------------------
% This code can be used to generate Figure 5.9
%----------------------------------------------------------------
clear all
close all
 
%________________Implementation OF FT Window/Sinc________________
M = 500;
t = (-M:M)*1e-3/5; 
E(450:550) = 1;E(1001)=0;
T = t(550)-t(450); 
 
index=300:700;
 
%---Figure 5.9(a)-----------------------------------------------------
 
area(t(index)*1e3,E(index)); axis([min(t(index))*1e3 max(t(index))*1e3 0 1.15])
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
xlabel('Time (ms)'); ylabel('Amplitude');
colormap(gray);
 
%---Figure 5.9(b)-----------------------------------------------------
index = 430:570;
df = 1/(max(t)-min(t));
f = (-M:M)*df;
Ef = T*fftshift(fft(E))/length(450:550); 
figure;
area(f(index),abs(Ef(index)));
axis([min(f(index)) max(f(index)) 0 .023])
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
xlabel('Frequency (Hz)'); ylabel('Amplitude');%grid on; 
colormap(gray);
 
%________________Implementation OF DFT________________
clear all;
% TIME DOMAIN SIGNAL
t = (-10:9)*1e-3; N=length(t);
En(1:N) = 1; 
 
%---Figure 5.9(c)-----------------------------------------------------
figure;
stem(t*1e3,En,'k','LineWidth',3); axis([min(t)*1.2e3 max(t)*1.2e3 0 1.25])
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('Time [ms]'); ylabel('s[n]');%grid on; 
 
%-----------FREQ DOMAIN SIGNAL---
dt = t(2)-t(1);
BWt = max(t)-min(t)+dt;
df = 1/BWt;
f = (-10:9)*df;
Efn = BWt*fftshift(fft(En))/length(En); 
 
%---Figure 5.9(d)-----------------------------------------------------
figure;
stem(f,abs(Efn),'k','LineWidth',3);
axis([min(f) max(f) 0 1.15])
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
xlabel('Frequency [Hz]'); ylabel('S[k]');%grid on; 
colormap(gray);hold on 
 
%-----this part for the sinc template
clear En2;
En2(91:110) = En; 
En2(200) = 0;
Efn2 = BWt*fftshift(fft(En2))/length(En); 
f2 = min(f):df/10:(min(f)+df/10*199);
plot(f2,abs(Efn2),'k-.','LineWidth',1);
axis([min(f2) max(f2) 0 .023]); hold off
 
%---------ZERO PADDING ----------
%TIME DOMAIN
clear En_zero;
En_zero(20:39) = En; 
En_zero(60) = 0;
dt = 1e-3;
t2 = dt*(-30:29);
 
%---Figure 5.9(e)-----------------------------------------------------
figure;
stem(t2*1e3,En_zero,'k','LineWidth',3); axis([-dt*30e3 dt*29e3 0 1.25])
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
xlabel('Time [ms]'); ylabel('s[n]');%grid on; 
 
%FREQUENCY DOMAIN
Efn2_zero = BWt*fftshift(fft(En_zero))/length(En); 
f2 = min(f):df/3:(min(f)+df/3*59);
 
%---Figure 5.9(f)-----------------------------------------------------
figure;
plot(f2,abs(Efn2_zero),'k-.','LineWidth',1);hold on
stem(f2,abs(Efn2_zero),'k','LineWidth',3);
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('Frequency [Hz]'); ylabel('S[k]');%grid on; 
axis([min(f2) max(f2) 0 0.023]); hold off
