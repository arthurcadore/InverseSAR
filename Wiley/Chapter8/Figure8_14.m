%------------------------------------------------------------
% This code can be used to generate Figure 8.14 
%------------------------------------------------------------
% This file requires the following files to be present in the same
% directory:
%
% scat_field.mat 

clear all
close all
 
%---Load the Scattered Field --------------------------------
load scat_field 
 
% Npulse = 128;  % number of pulses in one burst          
% Nburst = 512;  % number of bursts 
% f1 = 3e9;      % starting frequency for the EM wave
% BWf = 512e6;   % bandwidth of the EM wave
% T1 = (Npulse-1)/BWf;      % pulse duration
% PRF = 20e3;    % Pulse Repetation Frequency
% PRI = 1/PRF;   % Pulse Repetation Interval
% W = 0.16;      % angular velocity  [rad/s] 
% Vr = 1.0;      % radial velocity  [m/s]
% ar = 0.0;      % acceleration [m/s2]
 
% c = 3.0e8;     % speed of the EM wave
 
%---Figure 8.14(a)--------------------------------------------
plot(-Xc,Yc,'square', 'MarkerSize',5,'MarkerFaceColor',[1 0 0]);
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
xlabel('X [m]'); 
ylabel('Y [m]');
axis([min(-Xc)*1.1 max(-Xc)*1.1 min(Yc)*1.1 max(Yc)*1.1])
 
%---Figure 8.14(b)-------------------------------------------
%---Form Classical ISAR Image -------------------------------
w=hanning(Npulse)*hanning(Nburst)';
Es=Es.*w;
Es(Npulse*4,Nburst*4)=0;
ISAR=abs(fftshift(ifft2((Es))));
figure;matplot2(XX,YY,ISAR,30);
colorbar; colormap(1-gray);
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
xlabel('Range [m]'); 
ylabel('X-Range [m]');
