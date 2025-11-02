%------------------------------------------------------------
% This code can be used to generate Figures 8.2 - 8.6
%------------------------------------------------------------
% This file requires the following files to be present in the same
% directory:
%
% Fighter.mat 

clear all
close all
clc

%---Radar parameters-----------------------------------------
pulses = 128;           % # no of pulses          
burst = 128;            % # no of bursts 
c = 3.0e8;              % speed of EM wave [m/s]
f0 = 10e9;      	      % Starting frequency of SFR radar system [Hz]
bw = 128e6;             % Frequency bandwidth [Hz]
T1 = (pulses-1)/bw;     % Pulse duration [s]
PRF = 20e3;             % Pulse repetition frequency [Hz]
T2 = 1/PRF;             % Pulse repetition interval [s]
dr = c/(2*bw);          % range resolution [m] 
 
%---Target parameters----------------------------------------
W = 0.03;   % Angular velocity [rad/s] 
Vr = 70.0;  % radial translational  velocity of EM wave [m/s]
ar = 0.1;   % radial accelation of EM wave [m/s^2]
R0 = 16e3;  % target's initial distance from radar [m]
theta0 = 0; % Initial angle of target's wrt target [degree]
 
%---Figure 8.2-----------------------------------------------
%load the coordinates of the scattering centers on the fighter
load Fighter 
 
h = plot(-Xc,Yc,'o', 'MarkerSize',8,'MarkerFaceColor',[1 0 0]);grid;
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
axis([-35 35 -30 30])
xlabel('X [m]'); ylabel('Y [m]');

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
        Es(burst,pulses) = 0.0;
        for scat = 1:1:length(Xc);     
            arg = (Tetw - theta(scat) );
            rngterm = R0 + Rvr - r(scat)*sin(arg);
            range = reshape(rngterm,pulses,burst);
            range = range.';
            phase = k_fac.* range;
            Ess = exp(j*phase);
            Es = Es+Ess;
        end
        Es = Es.';
 
%---Figure 8.3-----------------------------------------------
%Form ISAR Image (no compansation)
X = -dr*((pulses)/2-1):dr:dr*pulses/2;Y=X/2;
ISAR = abs(fftshift(fft2((Es))));
h = figure;
matplot2(X,1:pulses,ISAR,20);
colormap(1-gray); colorbar;
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
xlabel('Range [m]'); 
ylabel('Doppler index');
   
%--Cross-Correlation Algorithm Starts here-------------------
RP=(ifft(Es)).';% Form Range Profiles
 
for l=1:burst; % Cross-correlation between RPn & RPref
    cr(l,:) = abs(ifft(fft(abs(RP(1,:))).* fft(abs(conj(RP(l,:))))));
    pk(l) = find((max(cr(l,:))== cr(l,:)));%Find max. ind. (range shift) range)
end
 
Spk = smooth((0:pulses-1),pk,0.8,'rlowess');%smoothing the delays 
RangeShifts = dr*pk;% range shifts
SmRangeShifts = dr*Spk;% range shifts
 
RangeDif = SmRangeShifts(2:pulses)-SmRangeShifts(1:pulses-1);%range differences
RangeDifAv =  mean(RangeDif);% average range differences 
 
T_burst = T(pulses+1)-T(1); % time between the bursts
Vr_Dif = (-RangeDif/T_burst);% estimated radial velocity from each RP
Vr_av = (RangeDifAv /T_burst);% estimated radial velocity (average)
 
%---Figure 8.4-----------------------------------------------
h = figure;
plot(i,RangeShifts,'LineWidth',2);hold
plot(i,SmRangeShifts,'-.k.','MarkerSize',4);hold
axis tight
legend('RP shifts','Smoothed RP shifts');
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
xlabel('range profile index'); 

%---Figure 8.5-----------------------------------------------
h = figure;
subplot(211);plot(RangeDif,'LineWidth',2);
axis([1 burst -.75 -.25 ])
set(gca,'FontName', 'Arial', 'FontSize',10,'FontWeight', 'Bold'); 
xlabel('Range profile index'); 
ylabel('Range differences [m] ')
 
subplot(212);
plot(Vr_Dif,'LineWidth',2);
axis([1 burst Vr-5 Vr+5 ])
set(gca,'FontName', 'Arial', 'FontSize',10,'FontWeight', 'Bold'); 
xlabel('Range profile index'); 
ylabel('Radial speed [m/s] ')
text(15,74,['Actual Speed = ',num2str(Vr),' m/s ,  Est. average speed = ',num2str(-Vr_av),' m/s']);

% Compansating the phase  
f = (f0+df);% frequency vector
T = reshape(T,pulses,burst); %prepare time matrix
F = f.'*ones(1,burst); %prepare frequency matrix
Es_comp = Es.*exp((j*4*pi*F/c).*(Vr_av*T));%Phase of E-field is compansated
 
%---Figure 8.6-----------------------------------------------
win = hanning(pulses)*hanning(burst).'; %prepare window
ISAR = abs((fft2((Es_comp.*win)))); % form the image
ISAR2 = ISAR(:,28:128);   
ISAR2(:,102:128)=ISAR(:,1:27);
h = figure;
matplot2(Y,1:pulses,ISAR2,20); % motion compansated ISAR image
colormap(1-gray);colorbar;
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
xlabel('Range [m]'); ylabel('Doppler index');
title('Motion compansated ISAR image')
