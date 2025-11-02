%----------------------------------------------------------------
% This code can be used to generate Figure 7.9 thru 7.13
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
XX = -dx*M/2:dx/4:-dx*M/2+dx/4*(4*M-1);
YY = -dy*N/2:dy/4:-dy*N/2+dy/4*(4*N-1);
 
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
 
%________________GET THE DATA________________________________________
load Es60 % load E-scattered
load  planorteta60_2_xyout.mat % load plane outline

%________________MATCHING PURSUIT_____________________________________
collectedData = zeros(200,3); %initilize scattering center info
ES = Es;
Power1 = sum(sum(Es).^2); % initial power of the data
axisX = X(1):dx/4:X(32); 
axisY = Y(1):dy/4:Y(64);
cosPhi = cos(PHI);
sinPhi = sin(PHI);
 
for N = 1:250; % extract 250 scattering centers
    Amax = 0;
    p1Max = zeros(size(ES));
    for Xn = axisX
        for Yn = axisY
            p1 = exp(-j*2*K.'*(cosPhi.*Xn+sinPhi.*Yn));
            A = sum(sum(ES.*p1))/(size(ES,1)*size(ES,2));
            if A > Amax
                Amax = A;
                collectedData(N,1:3) = [A Xn Yn];
                p1Max = conj(p1);
            end
        end
    end
    ES = ES-(Amax.*p1Max);
end
 
%--------Field Reconsctruction-----------
Esr = zeros(32,64);
for N = 1:250
    A = collectedData(N,1);
    x1 = collectedData(N,2);
    y1 = collectedData(N,3);
    Esr = Esr+A*exp(j*2*K.'*(cosPhi.*x1+sinPhi.*y1));
end
 
%---Figure 7.9(a)------------------------------------------------
%---SCATTERING CENTER INFO DISPLAY-------------------------------
load  planorteta60_2_xyout.mat
SSs = collectedData;
h = figure; 
plot(abs(SSs(1:250,1)),'square','MarkerSize',4,'MarkerFaceColor',[1 0 0]);
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
xlabel('Scattering Center #'); 
ylabel('Amplitude [mV/m]'); 

%---Figure 7.9(b)------------------------------------------------
h = figure; 
hold
for m=1:150
    t = round(abs(SSs(m,1))*20/abs(SSs(1,1)))+1
   plot(-SSs(m,2),SSs(m,3),'o','MarkerSize',t,'MarkerFaceColor',[1 0 0]);
end
hold
line( -xyout_xout,xyout_yout,'Color','k','LineStyle','.');
axis([min(X) max(X) min(Y) max(Y)])
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
xlabel('Range [m]'); 
ylabel('X-Range [m]');
title('Locations of scattering centers with relative amplitudes ')
            
%--ISAR IMAGE COMPARISON------------------------------------------
Enew = Es;
Enew(M*4,N*4) = 0;
ISARorig = fftshift(fft2(Enew)); 
ISARorig = ISARorig/M/N;
 
h = figure; 
matplot(X,Y,abs(ISARorig(4*M:-1:1,:).'),20); 
colorbar;
colormap(1-gray); 
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
line(-xyout_xout,xyout_yout,'Color','k','LineStyle','.');
xlabel('Range [m]'); 
ylabel('X-Range [m]'); 
title('Original ISAR image')

Enew = Esr;
Enew(M*4,N*4) = 0;
ISARrec = fftshift(fft2(Enew)); 
ISARrec = ISARrec.'/M/N; % reconstructed ISAR image
 
%---Figure 7.13 --------------------------------------------------
h = figure; 
matplot(X,Y,abs(ISARrec(:,4*M:-1:1)),20); 
colorbar;
colormap(1-gray); 
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
line(-xyout_xout,xyout_yout,'Color','k','LineStyle','.');
xlabel('Range [m]'); 
ylabel('X-Range [m]'); 
title('Reconstructed ISAR image')

%-------FIELD COMPARISON---------------------
%-------------------------------------------
h = figure; 
matplot(F,PHI*180/pi,abs((Es.')),20); 
colorbar; 
colormap(1-gray)
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
ylabel('Angle [Degree]'); 
xlabel('Frequency [GHz]'); 
title('Original back-scattered field')

%---Figure 7.10 --------------------------------------------------
h = figure;
matplot(F,PHI*180/pi,abs((Esr.')),20); 
colorbar; 
colormap(1-gray)
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
ylabel('Angle [Degree]'); 
xlabel('Frequency [GHz]'); 
title('Reconstructed back-scattered field')

%--------RECONSTRUCT THE FIELD PATTERN x10-------------
Esr = zeros(320,640);
k = K;
kk = k(1):(k(32)-k(1))/319:k(32);
pp = PHI(1):(PHI(64)-PHI(1))/639:PHI(64);
csP = cos(pp); 
snP = sin(pp);
for N = 1:250
    A = collectedData(N,1);
    x1 = collectedData(N,2);
    y1 = collectedData(N,3);
    Esr = Esr+A*exp(j*2*kk.'*(csP.*x1+snP.*y1));
end
 
%---Figure 7.11 --------------------------------------------------
h = figure;
matplot(F,PHI*180/pi,abs((Esr.')),20); 
colorbar; 
colormap(1-gray)
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
ylabel('Angle [Degree]'); 
xlabel('Frequency [GHz]'); 
title('Reconstructed field (x10 upsampled)')

%---Figure 7.12(a)------------------------------------------------
nn = 7;
h = figure;
plot(PHI*180/pi,abs(Es(nn,:)),'k-.*','MarkerSize',8,'LineWidth',2);
hold;
plot(pp*180/pi,abs(Esr(10*(nn-1)+1,:)),'k-','LineWidth',2);
hold;
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
xlabel('PHI [Degree]'); ylabel('Scat. field [V/m]'); 
tt = num2str(F(nn)); 
ZZ = ['@ f =  ' tt '  GHz'];
axis([PHI(1)*180/pi PHI(64)*180/pi 0 0.35])
title(ZZ);
drawnow; 
legend('with brute force computation','with scattering centers')

%---Figure 7.12(b)------------------------------------------------
nn = 4;
h = figure;
plot(F,abs(Es(:,nn)),'k-.*','MarkerSize',8,'LineWidth',2);
hold;
plot(kk*c/2/pi,abs(Esr(:,10*(nn-1)+1)),'k-','LineWidth',2);
hold;
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Bold');
xlabel('Frequency [GHz]'); ylabel('Scat. field [V/m]'); 
tt = num2str(PHI(nn)); 
ZZ = ['@ PHI =  ' tt '  Deg.'];
title(ZZ);
drawnow; 
axis([F(1) F(32) 0 0.5])
legend('with brute force computation','with scattering centers')
