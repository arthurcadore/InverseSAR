%----------------------------------------------------------------
% This code can be used to generate Figure 5-11
%----------------------------------------------------------------
clear all
close all
clc

% Prepare mesh
[X,Y] = meshgrid(-6:.1:6, -6:.1:6); 
M = length(X); 
N = length(Y) ; 
Object=zeros(M,N);

% Set 3 scattering centers
hh = figure;
Object(101,95)=5; 
Object(30,96)=2; 
Object(100,15)=3;
%---Figure 5-11(a)-----------------------------------------------------
surf(X,Y,Object); 
colormap(1-gray);
axis tight; 
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('X [m]'); 
ylabel('Y [m]');
zlabel('Amplitude') 
view(-45,20)
saveas(hh,'Figure5-11a.bmp','bmp');

%Find spatial resolutions
% fc = 10; % center frequency
% phic = 0; % center angle
% c = .3; % speed of light
dx = X(1,2)-X(1,1); % range resolution
dy = dx; % xrange resolution

%Find Bandwidth in spatial frequencies
BWkx = 1/dx; 
BWky = 1/dy; 

% PSF
h = sinc(BWkx*X/pi).*sinc(BWky*Y/pi);

%---Figure 5-11(b)-----------------------------------------------------
hh = figure;
surf(X,Y,abs(h)); 
axis tight; 
colormap(1-gray);
axis([-6 6 -6 6 0 1])
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('X [m]'); 
ylabel('Y [m]');
zlabel('Amplitude'); 
view(-45,20)
saveas(hh,'Figure5-11b.bmp','bmp');

%Convolution
hh = figure;
ISAR = fft2(fft2(Object).*fft2(h))/M/N;
%---Figure 5-11(c)-----------------------------------------------------
surf(X,Y,abs(ISAR)); 
axis tight; 
colormap(1-gray);
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('Range [m]'); 
ylabel('Cross - range [m]');
zlabel('ISAR ');
view(-45,20)
saveas(hh,'Figure5-11c.bmp','bmp');

%---Figure 5-11(c)-----------------------------------------------------
hh = figure;
matplot(X(1,1:M),Y(1:N,1),ISAR,30); 
colormap(1-gray);
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('Range [m]'); 
ylabel('Cross - range [m]');
title('ISAR ');
saveas(hh,'Figure5-11d.bmp','bmp');

