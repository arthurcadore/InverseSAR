%----------------------------------------------------------------
% This code can be used to generate Figure 3-14
%----------------------------------------------------------------
clear all
close all

%--- transmitted signal ------------------------------------
fc = 8e8; % initial frequency 
BWf = 10e6; % frequency bandwidth
To = 5e-6; %pulse duration
Beta = BWf/To; 
N = 400; %sample points
td = 1e-6; % delay 

t = -To/2:To/(N-1):To/2; %time vector
tt = t*1e6; %time vector in micro seconds
f  =fc:BWf/(N-1):(fc+BWf);% frequency vector 

s = 10*cos(2*pi*(fc*(t-td)+Beta*((t-td).^2)));% transmitted signal
sr = 10*cos(2*pi*(fc*t+Beta*(t.^2)));% replica

%---Figure 3.14(a)------------------------------------------------
h=figure;plot(tt,s, 'k','LineWidth',2)
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
title('transmitted signal');
xlabel(' Time [\mus]')
axis([min(tt) max(tt) -20 20 ]);
%saveas(h,'MF_1.png','png');

%---Figure 3.14(b)------------------------------------------------
%--- Noise Signal -----------------------------------------------
n=5*randn(1,N);
% Plot noise signal 
h=figure;plot(tt,n, 'k','LineWidth',2)
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
xlabel(' Time [\mus]'),title('noise signal '); 
axis([min(tt) max(tt) -20 20 ]);          
%saveas(h,'MF_2.png','png');

%---Figure 3.14(c)------------------------------------------------
%--- Received Signal --------------------------------------------
x=s+n;                
% Plot received signal 
h=figure;plot(tt,x, 'k','LineWidth',2)
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel(' Time [\mus]'),title('received signal'); 
axis([min(tt) max(tt) -20 20 ]);
%saveas(h,'MF_3.png','png');

%---Figure 3.14(d)------------------------------------------------
%--- Matched Filtering ------------------------------------------
X=fft(x)/N; 
S = conj(fft(sr)/N);  
H=S;  
Y=X.*H; 
y=fftshift(ifft(Y)); 
%----Plot matched filter output---------------------------------- 
h=figure;plot(tt,abs(y), 'k','LineWidth',2)
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel(' Time [\mus]')
axis([min(tt) max(tt) 0 .2]);title('Matched filter output ');
%saveas(h,'MF_4.png','png');
