%----------------------------------------------------------------
% This code can be used to generate Figure 1_8
%----------------------------------------------------------------

clear all
close all

%% DEFINE PARAMETERS
t=linspace(-50,50,1001); % Form time vector
df=1/(t(2)-t(1)); %Find frequency resolution
f=df*linspace(-50,50,1001);% Form frequency vector

%% FORM AND PLOT RECTANGULAR WINDOW
b(350:650)=ones(1,301);
b(1001)=0;
subplot(221);h=area(t,b); 
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
xlabel('Time [s]'); 
axis([-50 50 0 1.25])
set(h,'FaceColor',[.5 .5 .5])

subplot(222); h=area(f,fftshift(abs(ifft(b))));
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
xlabel('Frequency [Hz]')
axis([-40 40 0 .4])
set(h,'FaceColor',[.5 .5 .5])
 
%% FORM AND PLOT HANNING WINDOW
bb=b;
bb(350:650)=hanning(301)';

subplot(223); h=area(t,bb);
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
xlabel('Time [s]'); 
axis([-50 50 0 1.25])
set(h,'FaceColor',[.5 .5 .5])

subplot(224); h=area(f,fftshift(abs(ifft(bb))));
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
xlabel('Frequency [Hz]')
axis([-40 40 0 .2])
set(h,'FaceColor',[.5 .5 .5])
