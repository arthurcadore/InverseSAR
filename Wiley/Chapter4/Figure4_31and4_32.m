%----------------------------------------------------------------
% This code can be used to generate Figure 4.31 and 4.32
%----------------------------------------------------------------
% This file requires the following files to be present in the same
% directory:
%
% E_field.mat 
% planorteta60_2_xyout.mat
% planorteta60xzout.mat
 
clear all
close all
 
c = .3; % speed of light
fc = 8; % center frequency
phic = 0*pi/180; % center of azimuth look angles
thc = 60*pi/180; % center of elevation look angles
 
%________________PRE PROCESSING OF ISAR________________
BWx = 12; % range extend
M = 64;  % range sampling
BWy = 16; % x-range1 extend
N = 128;  % x-range1 sampling
BWz = 6; % x-range2 extend
P = 32;  % x-range2 sampling
 
%Find spatial resolutions
dx = BWx/M; % range resolution
dy = BWy/N; % xrange1 resolution
dz = BWz/P; % xrange1 resolution
 
% Form spatial vectors
X = -dx*M/2:dx:dx*(M/2-1);% range vector
XX = -dx*M/2:dx/4:-dx*M/2+dx/4*(4*M-1); % range vector (4x upsampled) 
Y = -dy*N/2:dy:dy*(N/2-1); % x-range1 vector
YY = -dy*N/2:dy/4:-dy*N/2+dy/4*(4*N-1); % x-range1 vector (4x upsampled) 
Z = -dz*P/2:dz:dz*(P/2-1); % x-range2 vector
ZZ = -dz*P/2:dz/4:-dz*P/2+dz/4*(4*P-1);% x-range2 vector (4x upsampled)
 
%Find resoltions in freq and angle
df = c/(2*BWx); 
dk = 2*pi*df/c; 
kc = 2*pi*fc/c;
 
dphi = pi/(kc*BWy);
dth = pi/(kc*BWz);
 
%Form F and PHI vectors
F = fc+[-df*M/2:df:df*(M/2-1)];
PHI = phic+[-dphi*N/2:dphi:dphi*(N/2-1)];
TET = thc+[-dth*P/2:dth:dth*(P/2-1)];
 
%________________GET THE DATA____________________________
load E_field % load E-scattered
load planorteta60_2_xyout % load target outline
 
% ISAR
ISAR=fftshift(ifftn(E3d)); 
 
%--------ISAR(x,y) Slices--------------
A = max(max(max(ISAR)));
for m=1:P;
    EE=ISAR(:,:,m); 
    EE(1,1,1)=A;
    zp = num2str(Z(m)); 
    zpp = ['ISAR(x,y) @ z = ' zp 'm'];
%---Figure 4.31-----------------------------------------------------
    matplot2(-X,Y, EE.',35);colorbar; colormap(1-gray);
    set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
    xlabel('Range [m]'); ylabel('X-Range [m]');
    drawnow; 
    title(zpp);
    h = line(-xyout_xout,xyout_yout,'Color','k','LineStyle','.','MarkerSize',5);
    pause
end
 
%---Figure 4.32-----------------------------------------------------
%--------XY Projection--------------
figure;
EExy = zeros(M,N);
load planorteta60_2_xyout
 
for m = 1:P;
    EExy = EExy+ISAR(:,:,m); 
end
    matplot2(-X,Y,EExy.',20);
    colorbar; 
    colormap(1-gray);
    set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
    xlabel('Range [m]'); 
    ylabel('X-Range [m]');
h=line(-xyout_xout,xyout_yout,'Color','k','LineStyle','.','MarkerSize',5);
  
 %--------XZ Projection--------------
figure;
load planorteta60xzout.mat
for m = 1:M;
    for n = 1:P
     EExz(m,n) = sum(ISAR(m,:,n)); 
    end
end
    matplot2(-X,Z,EExz.',20);colorbar; 
    colormap(1-gray);
    set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
    xlabel('Range [m]'); 
    ylabel('X-Range [m]');
h=line(-xzout_xout,-xzout_zout,'Color','k','LineStyle','.','MarkerSize',5);
     
  %--------YZ Projection--------------
figure;
load planorteta60yzout.mat
for m = 1:N;
    for n = 1:P
     EEyz(m,n) = sum(ISAR(:,m,n)); 
    end
end
    matplot2(-Y,Z,EEyz.',20);colorbar; 
    colormap(1-gray);
    set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
    xlabel('X-Range [m]'); 
    ylabel('X-Range [m]');
h=line(-yzout_yout,-yzout_zout,'Color','k','LineStyle','.','MarkerSize',5);
