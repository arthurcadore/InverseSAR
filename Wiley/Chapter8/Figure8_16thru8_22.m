%------------------------------------------------------------
% This code can be used to generate Figures 8.16 thru 8.22
%------------------------------------------------------------
% This file requires the following files to be present in the same
% directory:
%
% Fighter3.mat 

clear all
close all
clc
%---Radar parameters-----------------------------------------
pulses = 128;           % # no of pulses          
burst = 512;            % # no of bursts 
c = 3.0e8;              % speed of EM wave [m/s]
f0 = 3e9;               % Starting frequency of SFR radar system [Hz]
bw = 384e6;             % Frequency bandwidth [Hz]
T1 = (pulses-1)/bw;     % Pulse duration [s]
PRF = 20e3;             % Pulse repetition frequency [Hz]
T2 = 1/PRF;             % Pulse repetition interval [s]
theta0 = 0;             % Initial angle of target's wrt target [degree]
W = 0.15;               % Angular velocity [rad/s] 
Vr = 35.0;              % radial translational  velocity of EM wave [m/s]
ar = -1.9;              % radial accelation of EM wave [m/s^2]
R0 = 1.3e3;             % target's initial distance from radar [m]
dr = c/(2*bw);          % range resolution [m]  
theta0 = -30;           % Look angle of the target
 
%---Figure 8.16 ---------------------------------------------
%load the coordinates of the scattering centers on the fighter
load Fighter3 
h = plot(-Xc,Yc,'o', 'MarkerSize',8,'MarkerFaceColor',[1 0 0]);
grid on;
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
axis([-20 20 -20 20])
xlabel('X [m]'); ylabel('Y [m]');
 
%Scattering centers in cylindirical coordinates
[theta,r] = cart2pol(Xc,Yc);
theta = theta+theta0*0.017455329; %add initial angle
 
i = 1:pulses*burst;
T = T1/2+2*R0/c+(i-1)*T2;%calculate time vector
Rvr = Vr*T+(0.5*ar)*(T.^2);%Range Displacement due to radial vel. & acc. 
Tetw = W*T;% Rotational Displacement due to angular vel. 
 
i = 1:pulses;
df = (i-1)*1/T1; % Frequency incrementation between pulses
k = (4*pi*(f0+df))/c; 
k_fac = ones(burst,1)*k; 
 
%Calculate backscattered E-field  
        Es(burst,pulses) = 0.0;
        for scat = 1:1:length(Xc);     
            arg = (Tetw - theta(scat) );
            rngterm = R0 + Rvr - r(scat)*sin(arg);
            range = reshape(rngterm,pulses,burst);
            range = range.';
            phase = k_fac.* range;
            Ess = exp(-j*phase);
            Es = Es+Ess;
        end
        Es = Es.';
        
%---Figure 8.17 ---------------------------------------------
%Form ISAR Image (no compansation)
X = -dr*((pulses)/2-1):dr:dr*pulses/2;Y=X/2;
ISAR = abs((fft2((Es))));
h = figure;
matplot2(X,1:pulses,ISAR,20);
colormap(1-gray); 
colorbar;
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
xlabel('Range [m]'); ylabel('Doppler index');
 
%---Figure 8.18 ---------------------------------------------
% JTF Representation of range cell  
EsMp = reshape(Es,1,pulses*burst);
S = spectrogram(EsMp,128,64,128);
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
 
% Searching the motion parameters via Matching Pursuit 
syc = 1;
RR = 1e3:1e2:2e3;
V = 10:40;
A = -2.5:.1:1;
m = 0; 
clear EI
 
for Vest = V;
    m = m+1;
    n = 0;
    for iv = A;
        n = n+1; 
        p = 0;
        for Rest = RR;
        p = p+1;
        VI(syc,1:2) = [Vest,iv];
        S = exp((j*4*pi*F/c).*(Rest+Vest*T+(0.5*iv)*(T.^2)));
        Scheck = Es.*S;
        SumU = sum(sum(Scheck));
        EI(m,n,p) = abs(SumU);
        end    
end
end
 
[dummy,pp] = max(max(max((EI)))); %Find index for estimated velocity
[dummy,nn] = max(max((EI(:,:,pp)))); %Find index for estimated velocity
[dummy,mm] = max(EI(:,nn,pp)); %Find index for estimated acceleration

%---Figure 8.19 ---------------------------------------------
figure;
h = surfc(A,V,EI(:,:,pp));
colormap(gray)
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
ylabel('Translational velocity [m/s]'); 
xlabel('Translational acceleration [m/s^2]');
zlabel ('maximum argument')
%saveas(h,'Figure9-16.png','png');
 
% Form the mathing phase for compensation
Sconj = exp((j*4*pi*F/c).*(V(mm)*T+(0.5*A(nn)*(T.^2))));
% Compansate
S_Duz = Es.*Sconj;
 
%---Figure 8.20 ---------------------------------------------
% ISAR after compensation 
h = figure;
matplot(X,burst,abs(fftshift(fft2(S_Duz))),30);
colormap(1-gray);
colorbar;%grid;
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
xlabel('Range [m]'); ylabel('Doppler index');

%---Figure 8.21 ---------------------------------------------
% Check the compensation using via JTF Representation of range cells 
Sconjres = reshape(S_Duz,1,pulses*burst);
S = spectrogram(Sconjres,128,64,120);
[a,b] = size(S);
h =figure; 
matplot2((1:a),(1:b),abs(S),50);
colormap(1-gray);
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
xlabel('time pulses'); ylabel('frequency index');
title('Spectrogram');
 
%---This part for Rotational motion compensation------------------
Ese = S_Duz;
win = hamming(pulses)* hamming(burst).';% Prepare Window
Esew = Ese.*win;                        % Apply window to the E-field
Es_IFFT = (ifft(Esew)).';               % Range profiles 
 
i = 1:pulses*burst;
T = T1/2+2*R0/c+(i-1)*T2;
 
%** Apply Gaussian Blur Filter via Gabor Function
N = 1;                        % Sampling #
tst = T2*pulses;              % dwell time
t = T(1:pulses:pulses*burst); % time vector for bursts
 
fp = 160;                   
Alpha_p = 0.04;              % Blurring coefficient
t_istenen = 100;             % tp=1 sec, T=2.1845 sec
tp = ((t_istenen-1)*tst)/2;  % Center of Gaussian window
 
% % Gabor function and Gaussian Blur function
    parca1 = 1/sqrt(2*pi*(Alpha_p)^2);  % normalized term
    parca2 = exp(-((t-tp).^2)/(2*Alpha_p)); % Gaussian window term
for i=1:pulses      
    parca3 = exp((-j*2*pi*fp*(t-tp))/N);  % Harmonic function
    GaborWavelet(i,1:burst) = parca1*parca2.*parca3;% Gabor Wavelet func.
    fp=fp+1/(pulses*tst);              
end;
 
% %** Wavelet Transform
St_Img = fftshift(GaborWavelet*Es_IFFT); % shift image to the center    
 
h = figure;
matplot2(X,pulses,(St_Img.'),30);
colorbar;
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
colormap(1-gray);
grid;
xlabel('Range [m]'); 
ylabel('Doppler index');
 
%---Figure 8.22 ---------------------------------------------
EMp = reshape(St_Img,1,128*128);
S = spectrogram(EMp,256,120);
h = figure;
matplot2((1:pulses),(1:pulses),abs(S.'),60);
colormap(1-gray);
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
xlabel('time pulses'); ylabel('frequency index');
title('Spectrogram');
