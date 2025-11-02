%----------------------------------------------------------------
% This code can be used to generate Figure 5.19 (a-b)
%----------------------------------------------------------------
% This file requires the following files to be present in the same
% directory:
%
% Esplanorteta60.mat 
% planorteta60_2_xyout.mat

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

% Load the backscattered data
load Esplanorteta60
load planorteta60_2_xyout

%________________ POST PROCESSING OF ISAR________________
ISAR = fftshift(fft2(Es.')); 
ISAR = ISAR/M/N;
%---Figure 5.19(c)-----------------------------------------------------
h = figure;
matplot2(X(32:-1:1),Y,ISAR,20); 
colormap(1-gray); 
colorbar
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('Range [m]'); 
ylabel('Cross - range [m]');%grid on; 
h = line(-xyout_xout,xyout_yout,'Color','k','LineStyle','.','MarkerSize',3);
%saveas(h,'Planor.png','png');

%windowing;
w = hamming(M)*hamming(N).';
Ess = Es.*w;
  
%zero padding;
Enew = Ess;
Enew(M*4,N*4) = 0;
XX = X(1):dx/4:X(1)+dx/4*(4*M-1);
YY = Y(1):dy/4:Y(1)+dy/4*(4*N-1);

% ISAR image formatiom
ISARnew = fftshift(fft2(Enew.')); 
ISARnew = ISARnew/M/N;

%---Figure 5.19(d)-----------------------------------------------------
load planorteta60_2_xyout.mat
h = figure;
matplot2(XX(4*M:-1:1),YY,abs(ISARnew),20);
colormap(1-gray); 
colorbar
line(-xyout_xout,xyout_yout,'Color','k','LineStyle','.','MarkerSize',3);
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('Range [m]'); 
ylabel('Cross - range [m]');
%saveas(h,'Planor4xHamming.png','png');

