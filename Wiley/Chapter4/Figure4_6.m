%----------------------------------------------------------------
% This code can be used to generate Figure 4-6
%----------------------------------------------------------------
clear all
close all

c = .3; % speed of light
fc = 4; % center frequency
phic = 0*pi/180; % center of azimuth look angles
thc = 80*pi/180; % center of elevation look angles

%________________PRE PROCESSING________________
BWx = 80; % range extend
M = 32;  % range sampling
dx = BWx/M; % range resolution
X = -dx*M/2:dx:dx*(M/2-1); % range vector
XX = -dx*M/2:dx/4:-dx*M/2+dx/4*(4*M-1); % range vector (4x upsampled) 

%Form frequency vector
df = c/2/BWx; % frequency resolution  
F = fc+[-df*M/2:df:df*(M/2-1)]; % frequency vector
k = 2*pi*F/c; % wavenumber vector

% load backscattered field data for the target
load Es_range

%zero padding (4x);
Enew=E;
Enew(M*4)=0;

% RANGE PROFILE GENERATION
RP = M*fftshift(ifft(Enew));
h = plot(XX,abs(RP),'k','LineWidth',2);
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
ylabel('range profile intensity'); xlabel('range [m]');
axis tight


 
