%----------------------------------------------------------------
% This code can be used to generate Figure 4-8
%----------------------------------------------------------------
clear all
close all

c = .3; % speed of light
fc = 4; % center frequency
phic = 0*pi/180; % center of azimuth look angles
thc = 80*pi/180; % center of elevation look angles

%________________PRE PROCESSING________________
BWy = 66; % x-range extend
N = 128;  % x-range sampling
dy = BWy/N; % x- range resolution
Y = -dy*N/2:dy:dy*(N/2-1); % x-range vector
YY = -dy*N/2:dy/4:-dy*N/2+dy/4*(4*N-1); % range vector (4x upsampled) 

%Form angle vector
kc = 2*pi*fc/c; % center wavenumber 
dphi = pi/(kc*BWy); % azimuth angle resolution
PHI=phic+[-dphi*N/2:dphi:dphi*(N/2-1)]; % azimuth angle vector

% load backscattered field data for the target
load Es_xrange

%zero padding (4x);
Enew = E;
Enew(N*4) = 0;

% X-RANGE PROFILE GENERATION
XRP = N*fftshift(ifft(Enew));
h = plot(YY,abs(XRP),'k','LineWidth',2);
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
ylabel('cross-range profile intensity'); xlabel('cross-range [m]');
axis tight


