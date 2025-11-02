%----------------------------------------------------------------
% This code can be used to generate Figure 3.10
%----------------------------------------------------------------
clear all
close all

T=1e-3;A=1;
tau=-2e-3:1e-5:2e-3;fd=-3e3:10:3e3;fd=fd.';
X1=sinc(fd*(T-abs(tau)));
TT=(T-abs(tau));
p=find(TT<0);
TT(p)=0;
X2=A*A*ones(length(fd),1)*TT;
X=X1.*X2; X=X/max(max(abs(X)));

%---Figure 3.10(a)------------------------------------------------
mesh(tau*1e3,fd*1e-3,abs(X));
colormap(gray)
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
xlabel ('Time Delay [ms]')
ylabel ('Doppler Shift [KHz]')
zlabel ('Normalized AF')

%---Figure 3.10(b)------------------------------------------------
imagesc(tau*1e3,fd*1e-3,abs(X));
colormap(gray)
colorbar
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight','Bold');
xlabel ('Time Delay [ms]')
ylabel ('Doppler Shift [KHz]')
