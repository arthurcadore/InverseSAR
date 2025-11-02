%----------------------------------------------------------------
% This code can be used to generate Figure 3.9
%----------------------------------------------------------------
clear all
close all

tau=-10:.1:10; fd=tau;L=length(fd);
dummy=ones(L,L);
ideal=fftshift(ifft2(dummy));
mesh(tau,fd,abs(ideal));
colormap(gray)
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
xlabel ('Time Delay')
ylabel ('Doppler Shift')
zlabel ('Normalized AF')

