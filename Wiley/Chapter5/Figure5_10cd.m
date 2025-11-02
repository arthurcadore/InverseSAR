%----------------------------------------------------------------
% This code can be used to generate Figure 5-10 (c-d)
%----------------------------------------------------------------
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
load Esairbus
load airbusteta80_2_xyout

%________________ POST PROCESSING OF ISAR________________
ISAR = fftshift(fft2(Es.'));
ISAR = ISAR/M/N;
%---Figure 5-10(c)-----------------------------------------------------
h=figure;
matplot2(X(32:-1:1),Y,ISAR,30); colormap(1-gray); colorbar
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
xlabel('Range [m]'); 
ylabel('Cross - range [m]');%grid on; 
colormap(1-gray);
line(-xyout_xout,xyout_yout,'Color','k','LineStyle','.');
saveas(h,'Figure5-10c.png','png');

%zero padding with 4 times;
Enew = Es;
Enew(M*4,N*4) = 0;

% ISAR image formatiom
ISARnew = fftshift(fft2(Enew.')); 
ISARnew = ISARnew/M/N;

%ISARnew(1,1)=2.62
load airbusteta80_2_xyout.mat;

%---Figure 5-10(d)-----------------------------------------------------
h=figure;
matplot2(X(32:-1:1),Y,ISARnew,30); 
colormap(1-gray); 
colorbar;
line(-xyout_xout,xyout_yout,'Color','k','LineStyle','.');
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
xlabel('Range [m]');
ylabel('Cross - range [m]');
saveas(h,'Figure5-10d.png','png');
