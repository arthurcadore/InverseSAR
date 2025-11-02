%----------------------------------------------------------------
% This code can be used to generate Figure 4.26 thru 4.28
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
%  WIDE BW AND WIDE ANGLE ISAR  
%------------------------------------------------
%  B- POLAR REFORMATTING
%------------------------------------------------
nSampling = 1500; % sampling number for integration
 
% Define Bandwidth
f = fMin:(fMax-fMin)/(nSampling):fMax;
k = 2*pi*f/.3;
kMax = max(k);
kMin = min(k);
 
% Define Angle
phi = phiMin:(phiMax-phiMin)/(nSampling):phiMax;
 
kc = (max(k)+min(k))/2;
 
kx=k.'*cos(phi);
ky=k.'*sin(phi);
 
kxMax = max(max(kx));
kxMin = min(min(kx));
kyMax = max(max(ky));
kyMin = min(min(ky));
 
MM=4; % up sampling ratio
clear kx ky;
kxSteps = (kxMax-kxMin)/(MM*(nSampling+1)-1);
kySteps = (kyMax-kyMin)/(MM*(nSampling+1)-1);
kx = kxMin:kxSteps:kxMax; Nx=length(kx);
ky = kyMin:kySteps:kyMax; Ny=length(ky);
kx(MM*(nSampling+1)+1) = 0;
ky(MM*(nSampling+1)+1) = 0;
 
%________________ FORM RAW BACKSCATTERED DATA________________
%load scattering centers
load fighterSC
l = length(xx);
 
%form backscattered E-field from scattering centers
Es = zeros((nSampling+1),(nSampling+1));
for n=1:length(xx);
    Es = Es+exp(-j*2*k.'*cos(phi)*xx(n)).*exp(-j*2*k.'*sin(phi)*yy(n));
end
 
%---Figure 4.24-----------------------------------------------------
matplot2(f,phi*180/pi,Es,40);
colormap(1-gray); colorbar
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('Frequency [GHz]'); 
ylabel('Angle [Degree]');
        
newEs = zeros(MM*(nSampling+1)+1,MM*(nSampling+1)+1);
t = 0;
v = 0;
for tmpk = k
    t = t+1;
    v = 0;
    for tmpPhi = phi
        v = v+1;
        tmpkx = tmpk*cos(tmpPhi);
        tmpky = tmpk*sin(tmpPhi);
        indexX = floor((tmpkx-kxMin)/kxSteps)+1;
        indexY = floor((tmpky-kyMin)/kySteps)+1;
        
        r1 = sqrt(abs(kx(indexX)-tmpkx)^2+abs(ky(indexY)-tmpky)^2);
        r2 = sqrt(abs(kx(indexX+1)-tmpkx)^2+abs(ky(indexY)-tmpky)^2);
        r3 = sqrt(abs(kx(indexX)-tmpkx)^2+abs(ky(indexY+1)-tmpky)^2);
        r4 = sqrt(abs(kx(indexX+1)-tmpkx)^2+abs(ky(indexY+1)-tmpky)^2);
        
        R = 1/r1+1/r2+1/r3+1/r4;
        
        A1 = Es(t,v)/(r1*R);        
        A2 = Es(t,v)/(r2*R);
        A3 = Es(t,v)/(r3*R);        
        A4 = Es(t,v)/(r4*R);
             newEs(indexY,indexX) = newEs(indexY,indexX)+A1;
             newEs(indexY,indexX+1) = newEs(indexY,indexX+1)+A2;
             newEs(indexY+1,indexX) = newEs(indexY+1,indexX)+A3;
             newEs(indexY+1,indexX+1) = newEs(indexY+1,indexX+1)+A4;
     end
end
 
% down sample newEs by MM times
newEs=newEs(1:MM: size(newEs),1:MM: size(newEs));
 
%---Figure 4.25-----------------------------------------------------
% reformatted data
h = figure;
Kx = kx(1:Nx-1); 
Ky = ky(1:Ny-1); 
matplot2(Kx,Ky,newEs,40); 
colormap(1-gray); 
colorbar
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('kx [rad/m]'); ylabel('ky [rad/m]');
 
% Find Corresponding ISAR window in Range and X-Range
kxMax = max(max(kx));
kxMin = min(min(kx));
kyMax = max(max(ky));
kyMin = min(min(ky));
 
BWKx = kxMax-kxMin;
BWKy = kyMax-kyMin;
 
dx = pi/BWKx; 
dy = pi/BWKy; 
X = dx*(-nSampling/2:nSampling/2);
Y = dy*(-nSampling/2:nSampling/2);
 
%---Figure 4.26-----------------------------------------------------
% Plot the resultant ISAR image
h = figure;
tt = nSampling/4:3*nSampling/4;
ISAR3 = fftshift(ifft2(newEs));
matplot2(X,Y,ISAR3(:,tt),25);  
axis([-8 8 -6 6])
colormap(1-gray);
colorbar
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('Range [m]'); 
ylabel('Cross - range [m]');
 
