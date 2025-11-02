%----------------------------------------------------------------
% This code can be used to generate Figure 6-21 thru23
%----------------------------------------------------------------
% This file requires the following files to be present in the same
% directory:
%
% CoutUssFletcher.mat 

clear all
close all
clc
 
%---Radar parameters------------------------------------------------
pulses = 128;             % # no of pulses          
burst = 128;              % # no of bursts 
c = 3.0e8;                % speed of EM wave [m/s]
f0 = 9e9;                 % Starting frequency of SFR radar system [Hz]
bw = 125e6;               % Frequency bandwidth [Hz]
T1 = (pulses-1)/bw;       % Pulse duration [s]
PRF = 35e3;               % Pulse repetition frequency [Hz]
T2 = 1/PRF;               % Pulse repetition interval [s]
 
%---target parameters------------------------------------------------
theta0 = 0;               % Initial angle of target's wrt target [degree]
w = 1.2;                  % Angular velocity [degree/s] 
Vr = 5.0;                 % radial velocity of EM wave [m/s]
ar = 0.04;                % radial accelation of EM wave [m/s^2]
R0 = 4e3;                 % target's initial distance from radar [m]
dr = c/(2*bw);            % range resolution [m]  
W = w*pi/180;             % Angular velocity [rad/s] 
 
%---load the coordinates of the scattering centers on the fighter------
load CoutUssFletcher 
 
%---Figure 6.21--------------------------------------------------------
n = 10; 
Xc =(xind(1:n:6142)-93.25)/1.2; 
Yc =-zind(1:n:6142)/1.2;
h = figure;
plot(Xc,-Yc,'o','MarkerSize',3,'MarkerFaceColor',[1 0 0]);
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
axis([-80 80 -60 90])
xlabel('X [m]'); 
ylabel('Y [m]');
 
%---Scattering centers in cylindirical coordinates------------------------
[theta,r] = cart2pol(Xc,Yc);
Theta = theta+theta0*0.017455329; %add initial angle
 
i = 1:pulses*burst;
T = T1/2+2*R0/c+(i-1)*T2;%calculate time vector
Rvr = Vr*T+(0.5*ar)*(T.^2);%Range Displacement due to radial vel. & acc. 
Tetw = W*T;% Rotational Displacement due to angular vel. 
      
i = 1:pulses;
df = (i-1)*1/T1; % Frequency incrementation between pulses
k = (4*pi*(f0+df))/c; 
k_fac = ones(burst,1)*k; 
 
%------Calculate backscattered E-field-------------------------------------
        Es(burst,pulses)=0.0;
        for scat=1:1:length(Xc);     
            arg = (Tetw - theta(scat) );
            rngterm = R0 + Rvr - r(scat)*sin(arg);
            range = reshape(rngterm,pulses,burst);
            range = range.';
            phase = k_fac.* range;
            Ess = exp(j*phase);
            Es = Es+Ess;
        end
        Es = Es.';
 
% define noise 
noise=10*randn(burst,pulses);
   
E_signal = sum(sum(abs(Es.^2)));
E_noise = sum(sum(abs(noise.^2)));
SNR = E_signal/E_noise
SNR_db = 10*log10(SNR)
 
Es = Es+noise.';
        
%---Figure 6.22---------------------------------------------------
% Check out the range profiles
X = -dr*((pulses)/2-1):dr:dr*pulses/2;Y=X/2;
RP = fft((Es.'));
RP = fftshift(RP,1);
h = figure;
matplot2(X,1:burst,RP.',20);
colormap(1-gray); 
colorbar;
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
xlabel('Range [m]'); 
ylabel('Burst index');
 
%Form ISAR Image (no compansation)
%---Figure 6.23---------------------------------------------------
ISAR = abs(fftshift(fft2((Es))));
h = figure;
matplot2(X,1:burst,ISAR(:,pulses:-1:1),25);
colormap(1-gray); 
colorbar;
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
xlabel('Range [m]'); 
ylabel('Doppler index');
