%----------------------------------------------------------------
% This code can be used to generate Figure 5-10 (a-b)
%----------------------------------------------------------------
clear all
close all

c=.3; % speed of light
%________________PRE PROCESSING OF ISAR________________
%Find spatial resolutions
BWx = 12;
BWy = 16;
M = 32; 
N = 64;
fc = 6;
phic=0;

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
k = 2*pi*F/c;

% Load the backscattered data
load Esplanorteta60
load planorteta60_2_xyout

%________________ POST PROCESSING OF ISAR________________
ISAR = fftshift(fft2(Es.')); 
ISAR = ISAR/M/N;

%---Figure 5-10(a)-----------------------------------------------------
h=figure;
matplot2(X(32:-1:1),Y,ISAR,20); colormap(1-gray); colorbar
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
xlabel('Range [m]'); ylabel('Cross - range [m]');%grid on; 
line(-xyout_xout,xyout_yout,'Color','k','LineStyle','.','MarkerSize',3);
saveas(h,'Figure5-10a.png','png');
%zero padding;
Enew = Es;
Enew(M*4,N*4) = 0;
XX = X(1):dx/4:X(1)+dx/4*(4*M-1);
YY = Y(1):dy/4:Y(1)+dy/4*(4*N-1);

% ISAR image formatiom
ISARnew = fftshift(fft2(Enew.')); 
ISARnew = ISARnew/M/N;

load planorteta60_2_xyout.mat
%---Figure 5-10(b)-----------------------------------------------------
h=figure;
matplot2(XX(4*M:-1:1),YY,abs(ISARnew),20);colormap(1-gray); colorbar
line(-xyout_xout,xyout_yout,'Color','k','LineStyle','.','MarkerSize',3);
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
xlabel('Range [m]'); ylabel('Cross - range [m]');%grid on; 
saveas(h,'Figure5-10b.png','png');

