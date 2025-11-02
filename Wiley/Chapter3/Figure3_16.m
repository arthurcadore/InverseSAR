%----------------------------------------------------------------
% This code can be used to generate Figure 3-16
%----------------------------------------------------------------
clear all
close all

Phi_3dB = 5 * pi/180 ; % 3dB beamwidth of the antenna : 5 degrees
R0=8e3; % radial distance of the scatterer 
f= 6e9; % frequency
lam=3e8/f; % wavelength

X_max = R0 *tan(Phi_3dB/2); % maximum cross-range extend
x=-X_max:2*X_max/99:X_max; % cross-range vector

R = R0 *(1+x.^2/R0^2).^(0.5); % real  range distance
R_est = R0+x.^2/2/R0; % estimated range distance

%---Figure 3.16(a)------------------------------------------------
h=figure;
plot(x,R/1e3,'k--','LineWidth',1); hold
plot(x,R_est/1e3,'k.','LineWidth',4);hold; grid on
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
legend('actual radial distance','estimated radial distance')
xlabel('synthetic Aperture [m]')
ylabel('distance [km]')
axis([min(x) max(x) R0/1e3-.25 R0/1e3+.25])

%---Figure 3.16(b)------------------------------------------------
h=figure;plot(x,(R-R_est)/lam,'k','LineWidth',2);grid
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
xlabel('synthetic Aperture [m]')
ylabel('range error value  [\lambda]')
axis([min(x) max(x) -1 1])

