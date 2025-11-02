%----------------------------------------------------------------
% This code can be used to generate Figure 4-15
%----------------------------------------------------------------
clear all
close all

c = .3; % speed of light
fc = 6; % center frequency

phic = 45*pi/180; % center of azimuth look angles

%________________PRE PROCESSING OF ISAR________________
BWx = 13; % range extend
M = 32; % range sampling
BWy = 13; % xrange extend 
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
K=2*pi*F/c; % wanenumber vector

%________________GET THE DATA____________________________
load PLANORPHI45_Es.mat; % load E-scattered
load planorphi45_2_xyout.mat; % load target outline

%________________ POST PROCESSING OF ISAR________________
%windowing;
w=hanning(M)*hanning(N).';
Ess=Es.*w;

%zero padding;
Enew=Ess;
Enew(M*4,N*4)=0;

% ISAR image formatiom
ISARnew=fftshift(ifft2(Enew)); 

h=figure;
matplot2(X,Y,abs(ISARnew),22); % form the image
colormap(1-gray);colorbar
line(xyout_yout,xyout_xout,'LineWidth',.25,'LineStyle','.','Color','k');
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('Range [m]'); ylabel('Cross-Range [m]'); 

