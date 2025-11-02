%----------------------------------------------------------------
% This code can be used to generate Figure 5.19 (e-f)
%----------------------------------------------------------------
% This file requires the following files to be present in the same
% directory:
%
% Esairbus.mat 
% airbusteta80_2_xyout.mat

clear all
close all

c=.3; % speed of light
%________________PRE PROCESSING OF ISAR________________
%Find spatial resolutions
BWx = 80; 
BWy = 66; 
M = 32; 
N = 64; 
fc = 4;
phic = 0;

% Image resolutions
dx = BWx/M; 
dy = BWy/N;

% Form spatial vectors
X = -dx*M/2:dx:dx*(M/2-1);
Y = -dy*N/2:dy:dy*(N/2-1);

%Find resoltions in freq and angle
df = c/(2*BWx); 
dk = 2*pi*df/c; 
kc = 2*pi*fc/c;
dphi = pi/(kc*BWy);

%Form F and PHI vectors
F = fc+[-df*M/2:df:df*(M/2-1)];
PHI = phic+[-dphi*N/2:dphi:dphi*(N/2-1)];

load Esairbus
load airbusteta80_2_xyout

%________________ POST PROCESSING OF ISAR________________
ISAR = fftshift(fft2(Es.')); 
ISAR = ISAR/M/N;
%---Figure 5.19(a)-----------------------------------------------------
h = figure;
matplot2(X(32:-1:1),Y,ISAR,30); 
colormap(1-gray); 
colorbar
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('Range [m]'); ylabel('Cross - range [m]');%grid on; 
colormap(1-gray);%colorbar
h=line(-xyout_xout,xyout_yout,'Color','k','LineStyle','.');
%saveas(h,'airbus4x.png','png');

%windowing;
w = hamming(M)*hamming(N).';
Ess = Es.*w;

%zero padding with 4 times;
Enew = Ess;
Enew(M*4,N*4) = 0;

% ISAR image formatiom
ISARnew = fftshift(fft2(Enew.')); 
ISARnew = ISARnew/M/N;

%---Figure 5.19(b)-----------------------------------------------------
h = figure;
matplot2(X(32:-1:1),Y,ISARnew,30); 
colormap(1-gray); 
colorbar
line(-xyout_xout,xyout_yout,'Color','k','LineStyle','.');
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('Range [m]'); 
ylabel('Cross - range [m]'); 
%saveas(h,'airbus4xHamming.png','png');