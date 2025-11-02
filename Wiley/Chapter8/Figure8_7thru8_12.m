%------------------------------------------------------------
% This code can be used to generate Figures 8.7 thru 8.12
%------------------------------------------------------------
% This file requires the following files to be present in the same
% directory:
%
% fighter2.mat 

clear all
close all
clc

%---Radar parameters-----------------------------------------
pulses = 128;           % # no of pulses          
burst = 128;            % # no of bursts 
c = 3.0e8;              % speed of EM wave [m/s]
f0 = 8e9;               % Starting frequency of SFR radar system [Hz]
bw = 384e6;             % Frequency bandwidth [Hz]
T1 = (pulses-1)/bw;     % Pulse duration [s]
PRF = 14.5e3;           % Pulse repetition frequency [Hz]
T2 = 1/PRF;             % Pulse repetition interval [s]
dr = c/(2*bw);          % slant range resolution [m] 
 
%---Target parameters----------------------------------------
W = 0.06;  % Angular velocity [rad/s] 
Vr = 4.0;  % radial translational  velocity of EM wave [m/s]
ar = 0.6;  % radial accelation of EM wave [m/s^2]
R0 =.5e3;  % target's initial distance from radar [m]
theta0 = 125; % Initial angle of target's wrt target [degree]
 
%---Figure 8.7-----------------------------------------------
%load the coordinates of the scattering centers on the fighter
load fighter2 
 
h = plot(Xc,Yc,'o', 'MarkerSize',8,'MarkerFaceColor',[1 0 0]);
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
axis([-20 20 -20 20])
xlabel('X [m]'); 
ylabel('Y [m]');

%Scattering centers in cylindirical coordinates
[theta,r]=cart2pol(Xc,Yc);
theta=theta+theta0*0.017455329; %add initial angle
 
i = 1:pulses*burst;
T = T1/2+2*R0/c+(i-1)*T2;%calculate time vector
Rvr = Vr*T+(0.5*ar)*(T.^2);%Range Displacement due to radial vel. & acc. 
Tetw = W*T;% Rotational Displacement due to angular vel. 
      
i = 1:pulses;
df = (i-1)*1/T1; % Frequency incrementation between pulses
k = (4*pi*(f0+df))/c; 
k_fac=ones(burst,1)*k; 
 
%Calculate backscattered E-field  
        Es(burst,pulses)=0.0; 
        for scat=1:1:length(Xc);     
            arg = (Tetw - theta(scat) );
            rngterm = R0 + Rvr - r(scat)*sin(arg);
            range = reshape(rngterm,pulses,burst);
            range = range.';
            phase = k_fac.* range;
            Ess = exp(-j*phase);
            Es = Es+Ess;
        end
        Es = Es.';
 
%---Figure 8.8 ----------------------------------------------
%Form ISAR Image (no compansation)
X = -dr*((pulses)/2-1):dr:dr*pulses/2;Y=X/2;
ISAR = abs(fftshift(fft2((Es))));
h = figure;
matplot2(X,1:pulses,ISAR,20);
colormap(1-gray); colorbar;
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
xlabel('Range [m]'); ylabel('Doppler index');
        
%---Figure 8.9 ----------------------------------------------
% JTF Representation of range cell  
EsMp = reshape(Es,1,pulses*burst);
S = spectrogram(EsMp,128,64,120);
[a,b] = size(S);

h = figure;
matplot2((1:a),(1:b),abs(S),50);
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
colormap(1-gray);
xlabel('time pulses');
ylabel('frequency index');
title('Spectrogram');

%Prepare time and frequency vectors
f = (f0+df);% frequency vector
T = reshape(T,pulses,burst); %prepare time matrix
F = f.'*ones(1,burst); %prepare frequency matrix
 
% Searching the motion parameters via min. entropy method  
syc=1;
V = -15:.2:15;
A = -0.4:.01:1;
m = 0; 
for Vest = V;
    m = m+1;
    n = 0;
    for iv = A;
        n = n+1;
        VI(syc,1:2) = [Vest,iv];
        S = exp((j*4*pi*F/c).*(Vest*T+(0.5*iv)*(T.^2)));
        Scheck = Es.*S;
        ISAR = abs(fftshift(fft2((Scheck))));
        SumU = sum(sum(ISAR));
        I = (ISAR/SumU);
        Emat = I.*log10(I);
        EI(m,n) = -(sum(sum(Emat)));
        syc = syc+1;
    end    
end
 
[dummy,mm] = min(min(EI.')); %Find index for estimated velocity
[dummy,nn] = min(min(EI));   %Find index for estimated acceleration
%---Figure 8_10 ---------------------------------------------
h =surfc(A,V,EI);
colormap(gray)
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
ylabel('Translational velocity [m/s]'); 
xlabel('Translational acceleration [m/s^2]');
zlabel ('Entropy value')
saveas(h,'Figure9-10.png','png');
 
% Form the mathing phase for compensation
Sconj = exp((j*4*pi*F/c).*(V(mm)*T+(0.5*A(nn)*(T.^2))));
% Compansate
S_Duz = Es.*Sconj;
 
%---Figure 8.11 ---------------------------------------------
% ISAR after compensation 
h = figure;
matplot2(X,burst,abs(fftshift(fft2(S_Duz))),20);
colormap(1-gray);
colorbar;%grid;
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
xlabel('Range [m]'); 
ylabel('Doppler index');

%---Figure 8.12 ---------------------------------------------
% Check the compensation using via JTF Representation of range cells 
EsMp = reshape(S_Duz,1,pulses*burst);
S = spectrogram(EsMp,128,64,120);
[a,b] = size(S);
h = figure;
matplot2((1:a),(1:b),abs(S),50);
colormap(1-gray);
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
xlabel('time pulses');
ylabel('frequency index');
title('Spectrogram');

