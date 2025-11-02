%----------------------------------------------------------------
% This code can be used to generate Figure 4-14
%----------------------------------------------------------------
clear all
close all

c = .3; % speed of light
fc = 10; % center frequency
phic = 180*pi/180; % center of azimuth look angles

%________________PRE PROCESSING OF ISAR________________
BWx = 3; % range extend
M = 16; % range sampling
BWy = 3; % xrange extend 
N = 32; % xrange sampling

dx = BWx/M; % range resolution
dy = BWy/N; % xrange resolution

% Form spatial vectors
X = -dx*M/2:dx:dx*(M/2-1);
Y = -dy*N/2:dy:dy*(N/2-1);

%Find resoltions in freq and angle
df = c/(2*BWx); % frequency resolution
dk = 2*pi*df/c; % wavenumber resolution
kc = 2*pi*fc/c;
dphi = pi/(kc*BWy);% azimuth resolution

%Form F and PHI vectors
F = fc+[-df*M/2:df:df*(M/2-1)]; % frequency vector
PHI = phic+[-dphi*N/2:dphi:dphi*(N/2-1)];% azimuth vector
K = 2*pi*F/c; % wavenumber vector

%________________GET THE DATA____________________________
load  Escorner

%________________POST PROCESSING OF ISAR________________

ISAR=fftshift(ifft2(Es)); 
h=figure;
matplot2(X,Y,abs(ISAR),20); % form the image
colormap (1-gray);
colorbar
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('Range [m]'); ylabel('Cross-Range [m]'); 
 
line([-0.7071 0], [-0.7071 0],'LineWidth',2,'Color','k');
line([-0.7071 0], [0.7071 0],'LineWidth',2,'Color','k');

