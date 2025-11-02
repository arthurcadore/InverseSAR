%----------------------------------------------------------------
% This code can be used to generate Figure 4.20
%----------------------------------------------------------------
% This file requires the following files to be present in the same
% directory:
%
% Esairbus.mat 
% airbusteta80_2_xyout.mat

clear all
close all
 
c = .3; % speed of light
fc = 4; % center frequency
phic = 0*pi/180; % center of azimuth look angles
thc = 80*pi/180; % center of elevation look angles
 
%________________PRE PROCESSING OF ISAR________________
BWx = 80; % range extend
M = 32;  % range sampling
BWy = 66; % x-range extend
N = 64;  % x-range sampling
 
 
dx = BWx/M; % range resolution
dy = BWy/N; % xrange resolution
 
% Form spatial vectors
X = -dx*M/2:dx:dx*(M/2-1);% range vector
XX = -dx*M/2:dx/4:-dx*M/2+dx/4*(4*M-1); % range vector (4x upsampled) 
Y = -dy*N/2:dy:dy*(N/2-1); % x-range vector
YY = -dy*N/2:dy/4:-dy*N/2+dy/4*(4*N-1); % range vector (4x upsampled) 
 
%________________GET THE DATA____________________________
load Esairbus   % load E-scattered
load airbusteta80_2_xyout.mat % load target outline
 
% ISAR 4x UPSAMPLED-------------------
%zero padding;
Enew = Es;
Enew(M*4,N*4)=0;
 
% ISAR image formatiom
h = figure;
ISARnew = fftshift(ifft2(Enew)); 
matplot2(X,Y,abs(ISARnew.'),30); % form the image
colormap(1-gray);colorbar
line(-xyout_xout,xyout_yout,'LineWidth',.25,'LineStyle','.','Color','k');
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
xlabel('Range [m]'); ylabel('Cross-Range [m]'); 
