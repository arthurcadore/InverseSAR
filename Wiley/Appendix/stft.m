function [B,T,F]=stft(Y,f,BW,r,d)
%STFT Calculates the Short Time Fourier Transform of Vector Y
 
% Inputs:
% Y  : signal in the frequency domain
% f  : vector of frequencies [Hz]
% BW : bandwidth (same unit as F) of the sliding window
% f  : desired dinamic range of the display
% d  : additional delay [s] (if desired)
%
% Outputs:
% B  : stft of vector Y
% F  : frequency vector [GHz]
% T  : frequency vector [ns]
 
% The window used is a Kaiser window with beta=6.0   
  
df = f(2)-f(1); % frequency resolution
Ws = round(BW/df); 
if (Ws<2),
  Ws=2;
end;
 
W = kaiser(Ws,6); % window is a Kaiser with beta=6.0
N = max(size(Y)); % find length of Y
 
% Spectrogram
[B,T,F] = specgram((Y.*exp(-j*2*pi*f*d))',N,1/df,W,Ws-1);
F = F+((Ws-1)/2)*df+f(1); % set frequency axis
T = T-d;  % set time axis
 
% Treshold the image to the dynamic range
bmax = max(max(abs(B)));
ra = bmax/(10^(r/20));
B = B.*(abs(B)>=ra)+ra*ones(size(B)).*(abs(B)<ra);
 
% Display the STFT
colormap(jet(256)); %set colormap
imagesc(T*10^(9),F*10^(-9),20*log10(abs(B')/bmax))
axis xy; % change origin location
xlabel('Time [ns]'), 
ylabel('Frequency [GHz]');
