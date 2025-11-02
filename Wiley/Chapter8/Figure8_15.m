%----------------------------------------------------------------
% This code can be used to generate Figure 8.15
%----------------------------------------------------------------
clear all
close all

%---Load the Scattered Field ------------------------------------
load scat_field 

% Npulse = 128;                       % number of pulses in one burst          
% Nburst = 512;                       % number of bursts 
% f1 = 3e9;                           % starting frequency for the EM wave
% BWf = 512e6;                        % bandwidth of the EM wave
% T1 = (Npulse-1)/BWf;                % pulse duration
% PRF = 20e3;                         % Pulse Repetation Frequency
% PRI = 1/PRF;                        % Pulse Repetation Interval
% W = 0.16;                           % angular velocity  [rad/s] 
% Vr = 1.0;                           % radial velocity  [m/s]
% ar = 0.0;                           % acceleration [m/s2]
c = 3.0e8;                            % speed of the EM wave
N=1;
T = 2*R0/c+((1:Npulse*Nburst)-1)*PRI; %
tst = PRI*Npulse;                     % Toplam Sinyal Tekrarlama Süresi
Nt=T(1:Npulse:Npulse*Nburst);

Es_IFFT = ifft(Es)';                  % take 1D IFFT of field  
%----Prepare JTF filter function------------------------ 
n=0;figure
for frame=90:115:1100;               % select time frames  
    n=n+1;                            % counter for plotting
    fp=145;                           % Görüntünün Ortalanamasý için Kullanýlan Frekans
    tp=((frame-1)*tst)/2              % window center for each frame
      
%----Prepare JTF filter function (Gabor function with Gaussian Blur)
    Alpha_p=(0.04);                   % Blurring coefficient
    for i=1:Npulse		
          part1=1/sqrt(2*pi*(Alpha_p)^2);                   % normalize terimi
          part2=exp(-((Nt-tp).^2)/(2*Alpha_p));              % Gaussian window terimi
          part3=exp((-j*2*pi*fp*(Nt-tp))/N);                 % Harmonik Fonksiyon
          GaborWavelet(i,1:Nburst)=part1*part2.*part3;     % Gaborun Elde Edilmesi
         fp=fp+1/(Npulse*tst);              
     end;
% %** Wavelet Transform
     St =fftshift(GaborWavelet*Es_IFFT);
     subplot(3,3,n);matplot2(XX,YY,(St.'),25); 
     colormap(1-gray);
     set(gca,'FontName', 'Arial', 'FontSize',10,'FontWeight','Bold');
     title(['t = ',num2str(tp),' s'],'FontAngle','Italic');
end
