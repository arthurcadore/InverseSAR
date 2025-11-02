%----------------------------------------------------------------
% This code can be used to generate Figure 3-8
%----------------------------------------------------------------
clear all
close all

%--- transmitted signal ------------------------------------
fc = 8e8; % initial frequency 
To = 10e-6; %pulse duration
N = 200; %sample points
td = 4e-6; % delay 

t = 0:To/(5*N-1):To; %time vector
tt = t*1e6; %time vector in micro seconds

s = 10*ones(1,N);
s(5*N) =0; % transmitted signal replica
sr = s;% replica
M=round(td/To*(5*N));% shift amount
ss=circshift(sr.',M);ss=ss.';

%---Figure 3.8(a)------------------------------------------------
h1=figure;
h = area(tt,sr);
set(h,'FaceColor',[.5 .5 .5])
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
title('transmitted signal');
xlabel(' Time [\mus]')
axis([min(tt) max(tt) -30 30]);
%saveas(h1,'MF_1.png','png');

%---Figure 3.8(b)------------------------------------------------
h1=figure;
h = area(tt,ss)
set(h,'FaceColor',[.5 .5 .5])
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
title('received signal without noise');
xlabel(' Time [\mus]')
axis([min(tt) max(tt) -30 30]);
%saveas(h1,'MF_2.png','png');

%---Figure 3.8(c)------------------------------------------------
%--- Noise Signal ------------------------------------
n=5*randn(1,5*N);
% Plot noise signal 
h1=figure;
h=area(tt,n)
set(h,'FaceColor',[.5 .5 .5]);
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
xlabel(' Time [\mus]'),title('noise signal '); 
axis([min(tt) max(tt) -30 30]);
%saveas(h1,'MF_3.png','png');

%---Figure 3.8(d)------------------------------------------------
%--- Received Signal ------------------------------------
x=ss+n;                
h1=figure;
h=area(tt,x)
set(h,'FaceColor',[.5 .5 .5]);
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel(' Time [\mus]'),title('received signal with noise'); 
axis([min(tt) max(tt) -30 30]);
%saveas(h1,'MF_4.png','png');

%---Figure 3.8(e)------------------------------------------------
%--- Matched Filtering ------------------------------------
X=fft(x)/N; 
S = conj(fft(sr)/N);  
H=S;  
Y=X.*H; 
y=ifft(Y); 

% Plot matched filter output 
h1=figure;
h=area(tt,real(y)); 
set(h,'FaceColor',[.5 .5 .5]);
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');xlabel(' Time [\mus]')
axis([min(tt) max(tt) -.1 2]); title('matched filter output '); 
%saveas(h1,'MF_5.png','png');
