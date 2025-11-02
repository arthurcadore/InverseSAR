%----------------------------------------------------------------
% This code can be used to generate Figure 4-21 and 4-22
%----------------------------------------------------------------
clear all
close all

c = .3; % speed of light
fc = 8; % center frequency
phic = 0*pi/180; % center of azimuth look angles

%________________PRE PROCESSING OF ISAR________________
BWx = 18; % range extend
M = 64; % range sampling
BWy = 16; % xrange extend 
N = 64; % xrange sampling

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

%________________ FORM RAW BACKSCATTERED DATA________________
%load scattering centers
load fighterSC
l = length(xx);
%---Figure 4.21-----------------------------------------------------
h = figure;
plot(xx,yy,'.')
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('Range [m]'); ylabel('Cross - range [m]');
colormap(1-gray);
xlabel('X [m]'); 
ylabel('Y [m]'); 
%saveas(h,'Figure4-21.png','png');

%form backscattered E-field from scattering centers
Es = zeros(M,N);
for m=1:l;
    Es = Es+1.0*exp(-j*2*K'*(cos(PHI)*xx(m)+sin(PHI)*yy(m)));
end

%_____ POST PROCESSING OF ISAR (small-BW small angles)________________
ISAR = fftshift(ifft2(Es.')); 

%---Figure 4.22-----------------------------------------------------
h = figure;
matplot2(X,Y,ISAR,25); colormap(1-gray); colorbar
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('Range [m]'); ylabel('Cross - range [m]');
colormap(1-gray);
%saveas(h,'Figure4-22.png','png');

