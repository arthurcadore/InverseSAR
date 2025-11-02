%----------------------------------------------------------------
% This code can be used to generate Figure 5.12 – 5.18
%----------------------------------------------------------------
% Comparison of windowing functions
%-----------------------------------
clear all
close all

N=33;

%---Figure 5.12(a)-----------------------------------------------------
%---Rectangular window
rect = rectwin(N);
h = figure;
area(rect);
grid;
colormap(gray)
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel(' samples '); 
ylabel('Amplitude'); 
title(' Rectangular Window')
axis([-2 N+2 0 2])
%saveas(h,'Rectangular.png','png');

%---Figure 5.12(b)-----------------------------------------------------
rect(16*N)=0;
Frect = fftshift(fft(rect)); 
Frect = Frect/max(abs(Frect));
h = figure;
plot(mag2db(abs(Frect)),'k','LineWidth',2);
grid
axis tight; 
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('samples '); 
ylabel('Normalized amplitude[dB]'); 
title ('Spectrum of Rectangular Window')
axis([1 16*N -120 3])
%saveas(h,'RectangularSpec.png','png');

%---Figure 5.13(a)-----------------------------------------------------
%---Triangular window
tri = triang(N);
h = figure;
area([0 tri.']);
grid; 
colormap(gray)
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel(' samples '); 
ylabel('Amplitude'); 
title (' Triangular Window')
axis([-2 N+4 0 2])
%saveas(h,'Triangular.png','png');

%---Figure 5.13(b)-----------------------------------------------------
tri(16*N)=0; 
Ftri = fftshift(fft(tri)); 
Ftri = Ftri/max(Ftri);
h = figure;
plot(mag2db(abs(Ftri)),'k','LineWidth',2); 
grid; 
hold off;
axis tight; set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('samples '); 
ylabel('Normalized amplitude [dB]'); 
title ('Spectrum of Triangular Window')
axis([1 16*N -120 3])
%saveas(h,'TriangularSpec.png','png');

%---Figure 5.14(a)-----------------------------------------------------
%---Hanning window
han = hanning(N);
h = figure;
area(han);
grid; 
colormap(gray);
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel(' samples '); 
ylabel('Amplitude'); 
title (' Hanning Window')
axis([-2 N+2 0 2])
%saveas(h,'Hanning.png','png');

%---Figure 5.14(b)-----------------------------------------------------
han(16*N) = 0;
Fhan = fftshift(fft(han)); 
Fhan = Fhan/max(Fhan);
h = figure;
plot(mag2db(abs(Fhan)),'k','LineWidth',2); 
grid; 
hold off;
axis tight; set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('samples ');
ylabel('Normalized amplitude [dB]'); 
title ('Spectrum of Hanning Window')
axis([1 16*N -120 3])
%saveas(h,'HanningSpec.png','png');

%---Figure 5.15(a)-----------------------------------------------------
%---Hamming window
ham = hamming(N); 
h = figure;
area(ham);
grid; 
colormap(gray);
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('samples '); 
ylabel('Amplitude'); 
title ('Hamming Window')
axis([-2 N+2 0 2])
%saveas(h,'Hamming.png','png');

%---Figure 5.15(b)-----------------------------------------------------
ham(16*N)=0;
Fham = fftshift(fft(ham)); 
Fham = Fham/max(Fham);
h = figure;
plot(mag2db(abs(Fham)),'k','LineWidth',2); 
grid; 
hold off;
axis tight; set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('samples '); 
ylabel('Normalized amplitude [dB]'); 
title ('Spectrum of Hamming Window')
axis([1 16*N -120 3])
%saveas(h,'HammingSpec.png','png');

%---Figure 5.16(a)-----------------------------------------------------
%---Kaiser window
ksr = kaiser(N,1.5*pi);
h = figure;
area(ksr);
grid; 
colormap(gray);
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('samples '); 
ylabel('Amplitude');
title ('Kaiser Window, Beta=1.5*pi')
axis([-2 N+2 0 2])
%saveas(h,'Kaiser.png','png');

%---Figure 5.16(b)-----------------------------------------------------
ksr(16*N) = 0; 
Fksr = fftshift(fft(ksr)); 
Fksr = Fksr/max(Fksr);
h = figure;
plot(mag2db(abs(Fksr)),'k','LineWidth',2); 
grid; 
hold off;
axis tight; set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('samples '); 
ylabel('Normalized amplitude [dB]');
title ('Spectrum of Kaiser Window, Beta=1.5*pi')
axis([1 16*N -120 3])
%saveas(h,'KaiserSpec.png','png');

%---Figure 5.17(a)-----------------------------------------------------
%---Blackman window
blk = blackman(N);
h = figure;
area(blk);
grid; 
colormap(gray);
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('samples '); 
ylabel('Amplitude'); 
title ('Blackman Window')
axis([-2 N+2 0 2])
%saveas(h,'Blackman.png','png');

%---Figure 5.17(b)-----------------------------------------------------
blk(16*N) = 0; 
Fblk = fftshift(fft(blk)); 
Fblk = Fblk/max(Fblk);
h = figure;
plot(mag2db(abs(Fblk)),'k','LineWidth',2);
grid; 
hold off;
axis tight; set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('samples ');
ylabel('Normalized amplitude [dB]'); 
title ('Spectrum of Blackman Window')
axis([1 16*N -120 3])
%saveas(h,'BlackmanSpec.png','png');

%---Figure 5.18(a)-----------------------------------------------------
%---Chebyshev window
cheby = chebwin(N);
h = figure;
area(blk);
grid; 
colormap(gray);
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('samples ');
ylabel('Amplitude'); 
title ('Chebyshev Window')
axis([-2 N+2 0 2])
% saveas(h,'Chebyshev.png','png');

%---Figure 5.18(b)-----------------------------------------------------
cheby(16*N)=0;
Fcheby = fftshift(fft(cheby)); 
Fcheby = Fcheby/max(Fcheby);
h = figure;
plot(mag2db(abs(Fcheby)),'k','LineWidth',2);
grid; 
hold off;
axis tight; set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
xlabel('samples '); 
ylabel('Normalized amplitude [dB]'); 
title ('Spectrum of Chebyshev Window')
axis([1 16*N -120 3])
%saveas(h,'ChebyshevSpec.png','png');


