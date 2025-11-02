%----------------------------------------------------------------
% This code can be used to generate Figure 4.23 and 4.24
%----------------------------------------------------------------
% This file requires the following files to be present in the same
% directory:
%
% fighterSC.mat 

clear all
close all
 
c = .3; % speed of light
fc = 8; % center frequency
fMin = 6; % lowest frequency
fMax = 10; % highest frequency
 
phic = 0*pi/180; % center of azimuth look angles
phiMin = -30*pi/180; % lowest angle
phiMax = 30*pi/180; % highest angle
%------------------------------------------------
%  WIDE-BW AND LARGE ANGLES ISAR 
%------------------------------------------------
%  A- INTEGRATION
%------------------------------------------------
nSampling = 300; % sampling number for integration
 
% Define Arrays
f = fMin:(fMax-fMin)/(nSampling-1):fMax;
k = 2*pi*f/.3; 
kMax = max(k); 
kMin = min(k);
kc = (max(k)+min(k))/2;
phi = phiMin:(phiMax-phiMin)/(nSampling-1):phiMax;
 
% resolutions
dx = pi/(max(k)-min(k)); % range resolution
dy = pi/kc/(max(phi*pi/180)-min(phi*pi/180)); % xrange resolution
 
% Form spatial vectors
X = -nSampling*dx/2:dx:nSampling*dx/2;
Y = -nSampling*dy/2:dy:nSampling*dy/2;
 
%________________ FORM RAW BACKSCATTERED DATA________________
%load scattering centers
load fighterSC
l = length(xx);
 
%form backscattered E-field from scattering centers
clear Es;
Es = zeros((nSampling),(nSampling));
for m=1:l;
    Es = Es+1.0*exp(-j*2*k.'*cos(phi)*xx(m)).*exp(-j*2*k.'*sin(phi)*yy(m));
end
 
axisX = min(xx)-1:0.05:max(xx)+1;
axisY = min(yy)-1:0.05:max(yy)+1;
 
% take a look at what happens when DFT is used 
%---Figure 4.23-----------------------------------------------------
ISAR1 = fftshift(ifft2(Es.')); 
matplot2(axisX,axisY,ISAR1,22); 
colormap(1-gray); colorbar
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('Range [m]'); ylabel('Cross - range [m]');
 
% INTEGRATION STARTS HERE
 
% Building Simpson Nodes; Sampling Rate is nSampling
% Weights over k
h = (kMax-kMin)/(nSampling-1);
k1 = (kMin:h:kMax).';
wk1 = ones(1,nSampling); 
wk1(2:2:nSampling-1) = 4; 
wk1(3:2:nSampling-2) = 2; 
wk1 = wk1*h/3;
 
% Weights over phi
h = (phiMax-phiMin)/(nSampling-1);
phi1 =  (phiMin:h:phiMax).';
wphi1 = ones(1,nSampling); 
wphi1(2:2:nSampling-1) = 4; 
wphi1(3:2:nSampling-2) = 2; 
 
wphi1 = wphi1*h/3;
% Combine for two dimensional integration
[phi1,k1] = meshgrid(phi1,k1); 
phi1 = phi1(:); 
k1 = k1(:);
w = wk1.'*wphi1; 
w = w(:).';
 
newEs = Es(:).';
newW = w.*newEs;
 
% Integrate
b = 2j;
ISAR2 = zeros((max(xx)-min(xx)+2)/0.05+1,(max(yy)-min(yy)+2)/0.05+1);
 
k1 = k1.*b;
cosPhi = cos(phi1);
sinPhi = sin(phi1);
 
tic;
x1 = 0;
for X1 = axisX
    x1 = x1+1;
    y1 = 0;
    for Y1 = axisY
        y1 = y1+1;
        ISAR2(x1,y1) = newW*(exp(k1.*(cosPhi.*X1+sinPhi.*Y1)));
    end
end
time1 = toc;
 
%---Figure 4.24-----------------------------------------------------
matplot2(axisX,axisY,ISAR2.',22); 
colormap(1-gray); colorbar
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('Range [m]'); 
ylabel('Cross - range [m]');
