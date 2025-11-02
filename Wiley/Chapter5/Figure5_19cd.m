%----------------------------------------------------------------
% This code can be used to generate Figure 5.19 (c-d)
%----------------------------------------------------------------
% This file requires the following files to be present in the same
% directory:
%
% ucak.mat 

clear all
close all

c=.3; % speed of light

%________________PRE PROCESSING OF ISAR________________
%Find spatial resolutions
BWx = 18; 
BWy = 16; 
M = 64; 
N = 64; 
fc = 8;
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
K = 2*pi*F/c;
%________________ FORM RAW BACKSCATTERED DATA________________
load ucak
l = length(xx);
Es = zeros(M,N);
for m=1:l;
    Es=Es+1.0*exp(j*2*K'*(cos(PHI)*xx(m)+sin(PHI)*yy(m)));
end

%_____ POST PROCESSING OF ISAR (Small BW Small angle)________________
ISAR = fftshift(fft2(Es.')); 
ISAR = ISAR/M/N;
h = figure;
matplot2(X(M:-1:1),Y,ISAR,25); 
colormap(1-gray); 
colorbar
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('Range [m]'); 
ylabel('Cross - range [m]');%grid on; 
colormap(1-gray);
% saveas(h,'Ucak.png','png');

%windowing;
w = hamming(M)*hamming(N).';
Ess = Es.*w;

%-------------zero padding with 4 times----------
Enew = Ess;
Enew(M*4,N*4) = 0;

% ISAR image formatiom
ISARnew = fftshift(fft2(Enew.')); 
ISARnew = ISARnew/M/N;

h = figure;
matplot2(X(M:-1:1),Y,ISARnew,25); 
colormap(1-gray); 
colorbar
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('Range [m]'); 
ylabel('Cross - range [m]');%grid on; 
% saveas(h,'Ucak4xHamming.png','png');

