%----------------------------------------------------------------
% This code can be used to generate Figure 7.2 thru 7.8
%----------------------------------------------------------------
% This file requires the following files to be present in the same
% directory:
%
% Es60.mat 
% planorteta60_2_xyout.mat

clear all
close all
 
c = .3; % speed of light
fc = 6; % center frequency
phic = 0*pi/180; % center of azimuth look angles
thc = 90*pi/180; % center of elevation look angles
 
%________________PRE PROCESSING OF ISAR________________
BWx = 12; % range extend
M = 32; % range sampling
BWy = 16; % xrange extend 
N = 64; % xrange sampling
 
dx = BWx/M; % range resolution
dy = BWy/N; % xrange resolution
 
% Form spatial vectors
X = -dx*M/2:dx:dx*(M/2-1);
Y = -dy*N/2:dy:dy*(N/2-1);
XX =-dx*M/2:dx/4:-dx*M/2+dx/4*(4*M-1);
YY =-dy*N/2:dy/4:-dy*N/2+dy/4*(4*N-1);
 
%Find resoltions in freq and angle
df = c/(2*BWx); % frequency resolution
dk = 2*pi*df/c; % wavenumber resolution
kc = 2*pi*fc/c;
dphi = pi/(kc*BWy);% azimuth resolution
 
%Form F and PHI vectors
F = fc+[-df*M/2:df:df*(M/2-1)]; % frequency vector
PHI = phic+[-dphi*N/2:dphi:dphi*(N/2-1)];% azimuth vector
K = 2*pi*F/c; % wavenumber vector
dk = K(2)-K(1); % wavenumber resolution
 
%________________GET THE DATA____________________________
load Es60 % load E-scattered
load  planorteta60_2_xyout.mat % load plane outline

% ISAR
ISAR = fftshift(fft2(Es)); 
ISAR = ISAR/M/N; % the image
 
% ISAR 4x UPSAMPLED-------------------
Enew = Es;
Enew(M*4,N*4) = 0;
ISARnew = fftshift(fft2(Enew)); 
ISARnew = ISARnew/M/N;
 
%________________2D-CLEAN____________________________
% prepare 2D sinc functions
 sincx = ones(1,M); 
 sincx(1,M+1:M*4) = 0; 
 hsncF = fft(sincx)/M;
 sincy=ones(1,N); 
 sincy(1,N+1:N*4) = 0; 
 hsncPHI = fft(sincy)/N;
 
 %initilize
 hh = zeros(4*M,4*N);
 
 ISARres = ISARnew.';
 Amax = max(max(ISARnew));
 ISARbuilt = zeros(N*4,M*4);
 
 % loop for CLEAN 
 for nn=1:250,
  [A,ix] = max(max(ISARres));
  [dum,iy] = max(max(ISARres.'));
  hsincX = shft(hsncF,ix);
  hsincY = shft(hsncPHI,iy);  
  hhsinc = hsincX.'*hsincY;
  ISARres = ISARres-A*hhsinc.';
  SSs(nn,1:3) = [A XX(ix) YY(iy)];
  II=ISARres; 
  II(1,1) = Amax;
  % Image Reconstruction
  ISARbuilt = ISARbuilt-A*hhsinc.';
 end
 
%________________IMAGE COMPARISON____________________________
%---Figure 7.2(a)---------------------------------------------------
h = figure; 
matplot(X,Y,abs(ISARnew(4*M:-1:1,:).'),20); 
colorbar;
colormap(1-gray); 
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
line( -xyout_xout,xyout_yout,'Color','k','LineStyle','.');
xlabel('Range [m]'); 
ylabel('X-Range [m]'); 
title('Original ISAR image')

%---Figure 7.4---------------------------------------------------
h = figure; 
matplot(X,Y,abs(ISARbuilt(:,4*M:-1:1)),20); 
colorbar;
colormap(1-gray); 
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
line( -xyout_xout,xyout_yout,'Color','k','LineStyle','.');
xlabel('Range [m]'); 
ylabel('X-Range [m]'); 
title('Reconstructed ISAR image')

%________________SCATTERING CENTER INFO DISPLAY_______________________
 
%---Figure 7.3-----------------------------------------------------
h = figure; 
plot(abs(SSs(:,1)),'square', 'MarkerSize',4,'MarkerFaceColor',[1 0 0]);
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
xlabel('Scattering Center #'); 
ylabel('Amplitude [mV/m]'); 

%---Figure 7.2(b)---------------------------------------------------
h = figure; 
hold
for m = 1:150
    t = round(abs(SSs(m,1))*20/abs(SSs(1,1)))+1
  plot(-SSs(m,2),SSs(m,3),'o', 'MarkerSize',t,'MarkerFaceColor',[1 0 0]);
end
hold
line(-xyout_xout,xyout_yout,'Color','k','LineStyle','.');
axis([min(X) max(X) min(Y) max(Y)])
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
xlabel('Range [m]'); 
ylabel('X-Range [m]');
title('Locations of scattering centers with relative amplitudes ')

 
%________________RECONSTRUCT THE FIELD PATTERN________________________
ESR = zeros(320,640);
ESr = zeros(32,64);
k = K;
kk = k(1):(k(32)-k(1))/319:k(32);
pp = PHI(1):(PHI(64)-PHI(1))/639:PHI(64);
for nn = 1:250;
    An = SSs(nn,1);    
    xn = SSs(nn,2);    
    yn = SSs(nn,3);
    ESR = ESR+An*exp(j*2*xn*(kk-k(1)).')*exp(j*2*kc*yn*(pp-PHI(1)));
    ESr = ESr+An*exp(j*2*xn*(k-k(1)).')*exp(j*2*kc*yn*(PHI-PHI(1)));
end
 
%---Figure 7.5(a)---------------------------------------------------
h = figure; 
matplot(F,PHI*180/pi,abs((Es.')),20); 
colorbar; 
colormap(1-gray)
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
ylabel('Angle [Degree]'); 
xlabel('Frequency [GHz]'); 
title('Original back-scattered field')

%---Figure 7.5(b)---------------------------------------------------
h = figure;
matplot(F,PHI*180/pi,abs((ESr.')),20); 
colorbar; 
colormap(1-gray)
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
ylabel('Angle [Degree]'); 
xlabel('Frequency [GHz]'); 
title('Reconstructed back-scattered field')

%---Figure 7.6-----------------------------------------------------
h = figure;
matplot(F,PHI*180/pi,abs((ESR.')),20); 
colorbar; 
colormap(1-gray)
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
ylabel('Angle [Degree]'); 
xlabel('Frequency [GHz]'); 
title('Reconstructed back-scattered field (x10 upsampled)')

%---Figure 7.7-----------------------------------------------------
nn = 3; 
h = figure;
plot(PHI*180/pi,abs(Es(nn,:)),'k-.*','MarkerSize',8,'LineWidth',2);
hold;
plot(pp*180/pi,abs(ESR(10*(nn-1)+1,:)),'k-','LineWidth',2);
hold;
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
xlabel('PHI [Degree]'); 
ylabel('Scat. field [V/m]'); 
tt=num2str(F(nn)); 
ZZ=['@ f =  ' tt '  GHz'];
axis([PHI(1)*180/pi PHI(64)*180/pi 0 0.5]);
title(ZZ);
drawnow; 
legend('with brute force computation','with scattering centers')

%---Figure 7.8-----------------------------------------------------
nn = 11;
h = figure;
plot(F,abs(Es(:,nn)),'k-.*','MarkerSize',8,'LineWidth',2);
hold;
plot(kk*c/2/pi,abs(ESR(:,10*(nn-1)+1)),'k-','LineWidth',2);
hold;
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
xlabel('Frequency [GHz]'); 
ylabel('Scat. field [V/m]'); 
tt=num2str(PHI(nn)); 
ZZ=['@ PHI =  ' tt '  Deg.'];
title(ZZ);drawnow; axis([F(1) F(32) 0 0.35])
legend('with brute force computation','with scattering centers')
