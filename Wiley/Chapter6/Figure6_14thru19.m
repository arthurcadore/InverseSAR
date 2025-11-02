%----------------------------------------------------------------
% This code can be used to generate Figure 6.14 thru 19
%----------------------------------------------------------------
% This file requires the following files to be present in the same
% directory:
%
% fighter.mat 

clear all
close all
 
%---Radar parameters------------------------------------------------
c = 3e8;                % speed of EM wave [m/s]
fc = 10e9;              % Center frequency of chirp [Hz]
BWf = 2.5e9;            % Frequency bandwidth of chirp [Hz] 
T1 =.4e-6;              % Pulse duration of single chirp [s]
 
%---target parameters------------------------------------------------
Vx = 120;               % radial translational  velocity of target [m/s]
Xo = 0e3;               % target's initial x coordinate wrt radar
Yo = 24e3;              % target's initial y coordinate wrt radar
Xsize = 180;            % target size in cross-range [m]
Ysize = 60;             % target size in range [m]
 
%----set parameters for ISAR imaging-------------------------------
%  range processing
Ro = sqrt(Xo^2+Yo^2);   % starting range distance [m]
dr = c/(2*BWf);         % range resolution [m]
fs = 2*BWf;             % sampling frequency [Hz]
M = round(T1*fs);       % range samples 
Rmax = M*dr;            % max. range extend [m]
RR = -M/2*dr:dr:dr*(M/2-1);     % range vector [m]
Xmax = 1*Xsize;         % range window in ISAR [m]
Ymax = 1*Ysize;         % cross-range window in ISAR [m]
 
% Chirp  processing
U = Vx/Ro;              % rotational velocity [rad/s]
BWdop = 2*U*Ysize*fc/c; % target Doppler bandwith [Hz]
PRFmin = BWdop;         % min. PRF
PRFmax = c/(2*Ro);      % max. PRF
N = floor(PRFmax/BWdop);% # of pulses
PRF = N*BWdop;          % Pulse repetition frequency [Hz]
T2 = 1/PRF;             % Pulse repetition interval 
T = T2*N;               % Dwell time [s] (also = N/PRF)
 
%  cross- range processing
dfdop = BWdop/N;        % doppler resolution
lmdc = c/fc;            % wavelength at fc 
drc = lmdc*dfdop/2/U;   % cross-range resolution
RC = -N/2*drc:drc:(N/2-1)*drc; % cross-range vector 
 
%---load the coordinates of the scattering centers on the fighter------
load fighter

%---Figure 6.14--------------------------------------------------------
h = figure;
plot(-Xc,Yc,'o', 'MarkerSize',8,'MarkerFaceColor',[1 0 0]);grid;
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
axis([-35 35 -30 30])
xlabel('X [m]'); ylabel('Y [m]');

%--- sampling & time parameters -----------------------------
dt = 1/fs;                  % sampling time interval 
t = -M/2*dt:dt:dt*(M/2-1);  % time vector along chirp pulse
XX = -Xmax/2:Xmax/(M-1):Xmax/2;
F = -fs/2:fs/(length(t)-1):fs/2;  % frequency vector
 
slow_t = -M/2*dt:T2: -M/2*dt+(N-1)*T2;
 
%--- transmitted signal ------------------------------------
Kchirp = BWf/T1;              % chirp pulse parameter
s = exp(j*2*pi*(fc*t+Kchirp/2*(t.^2)));% original signal
sr =exp(j*2*pi*(fc*t+Kchirp/2*(t.^2)));% replica
H = conj(fft(sr)/M);   % matched filter transfer function
 
%---Figure 6.15--------------------------------------------------------
h = figure;
plot(t*1e6, s, 'k','LineWidth',0.5)
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
title('transmitted signal');
xlabel(' Time [\mus]')
axis([min(t)*1e6 max(t)*1e6 -2 2 ]);
 
%---Figure 6.16--------------------------------------------------------
h = figure;plot(F*1e-9,abs(fftshift(H)), 'k','LineWidth',2)
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
title('Matched filter response');
xlabel(' Frequency [GHz]')
 
%--- Received Signal ------------------------------------
for n=1: N
    Es(n,1:M) =zeros(1,M);
    for m =1: length(Xc);
        x = Xo+Xc(m)-Vx*T2*(n-1);
        R = sqrt((Yo+Yc(m))^2+x^2); 
        Es(n,1:M) = Es(n,1:M)+exp(j*2*pi*(fc*(t-2*R/c)+Kchirp/2*((t-2*R/c).^2)));
    end
       % define noise 
       noise=5*randn(1,M);
       NS(n,1:M)=noise;
       
       % Matched filtering
    EsF(n,1:M) = fft(Es(n,1:M)+noise)/M;
    ESS(n,1:M) = EsF(n,1:M).*H;
    ESS(n,1:M) = ifft(ESS(n,1:M));
end;
   
E_signal = sum(sum(abs(EsF.^2)));
E_noise = sum(sum(abs(NS.^2)));
SNR = E_signal/E_noise;
SNR_db = 10*log10(SNR);
 
%---Figure 6.17--------------------------------------------------------
rd=30; % dynamic range of display
h=figure; 
matplot2(1:N,slow_t,(ESS),rd);
colormap(1-gray);
colorbar
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
xlabel('Range bins'); 
ylabel('Azimuth time [s]'); 
 
%---Figure 6.18--------------------------------------------------------
h=figure; 
matplot2(RC,1:N,fftshift(fft(ESS.*win)),rd);
colormap(1-gray);
colorbar
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
xlabel('Range [m]'); 
ylabel('Doppler index');
 
%---Figure 6.19--------------------------------------------------------
win = hanning(N)*ones(1,M); % prepare window in cross-range direction
h=figure; 
matplot2(RC,XX,fftshift(fft(ESS.*win)),rd);
colormap(1-gray);
colorbar
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
xlabel('Range [m]'); 
ylabel('Cross-Range [m]');
axis([-30 30 -60 60])

