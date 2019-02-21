clear all;  close all;  clc;

display('Running test');

xS    = [-100,0];  % source position: 100 m above center 
xR    = [-60,0];   % central point of receiver array (-50*dxR:50*dxR)
c_0   = 1500;      % wave speed in embedding
c_s   = 3000;      % wave speed in scattering object
rho_0 = 3000;      % mass density of enbedding
rho_s = 1500;      % mass density of scattering object

dxR   = 0.5;       % receiver grid size
gsize = [1+240/dxR,1+200/dxR]; % gridsize in [z,x], center of model at coordinate (0,0)
f_in  = 50;        % selected frequency to compare, returned in Pf
ntrcv = 256;       % number of time samples in FD modeling
dtrcv = 0.004;     % dt in receivers


% Make grid
a  = 40;                               % radius circle cylinder 
N1 = gsize(1);                         % number of samples in x_1  
N2 = gsize(2);                         % number of samples in x_2 
dx = dxR;                              % with meshsize dx
x1 = -(N1+1)*dx/2 + (1:N1)*dx;
x2 = -(N2+1)*dx/2 + (1:N2)*dx;
[X1,X2] = ndgrid(x1,x2);
% Now array subscripts are equivalent with Cartesian coordinates
% x1 axis points downwards and x2 axis is in horizontal direction
% x1 = X1(:,1) is a column vector in vertical direction
% x2 = X2(1,:) is a row vector in horizontal direction
R = sqrt(X1.^2 + X2.^2);
cgrid  = c_s * (R < a) + c_0 * (R >= a);  
rhogrid  = rho_s * (R < a) + rho_0 * (R >= a); 
   
% DATA from Thorbecke's finite difference code
[Ptsct, Pfinc, Pfsct, f_out]=FD_mod_grid( xS, xR, ntrcv, dtrcv, dxR, cgrid, rhogrid, f_in );

f=f_out           % nearest computed discrete frequency

% Compare with analytical solution ---------------------------------------
[P_inc,P_sct,xR] = ForwardCircle( xS, xR, gsize, dxR, c_0, c_s, rho_0, rho_s, f );

set(figure(1),'Units','centimeters','Position',[1 1 30 10]);
subplot(1,3,1)
  plot(xR(2,:),real(Pfinc),'LineWidth',1.2); 
  title('Real part of P^{inc}'); axis tight; hold on;
  plot(xR(2,:),real(P_inc),'--r','LineWidth',1.2); 
  axis tight; hold off;
subplot(1,3,2)
  plot(xR(2,:),imag(Pfinc),'LineWidth',1.2); 
  title('Imaginary part of P^{inc}'); axis tight; hold on;
  plot(xR(2,:),imag(P_inc),'--r','LineWidth',1.2);
  axis tight; hold off;
subplot(1,3,3)
  plot(xR(2,:), abs(Pfinc),'LineWidth',1.2); 
  title('Absolute value of P^{inc}'); axis tight; hold on; 
  plot(xR(2,:), abs(P_inc),'--r','LineWidth',1.2); 
  axis tight; hold off
legendtitle1 = sprintf('f=%.2fHz FiniteDiff', f);
legendtitle2 = sprintf('f=%.2fHz Analytic  ', f);
legend(legendtitle1,legendtitle2,'Location','Best');

set(figure(2),'Units','centimeters','Position',[9 12 10 10]);
error = abs(P_sct(:)-Pfsct(:))./abs(Pfsct(:));
plot(error,'LineWidth',1.2); title('Relative error'); axis tight;

set(figure(3),'Units','centimeters','Position',[1 1 30 10]);
subplot(1,3,1)
  plot(xR(2,:),real(Pfsct),'LineWidth',1.2); 
  title('Real part of P_{sct}'); axis tight; hold on;
  plot(xR(2,:),real(P_sct),'LineWidth',1.2); 
subplot(1,3,2)
  plot(xR(2,:),imag(Pfsct),'LineWidth',1.2); 
  title('Imaginary part of P^{sct}'); axis tight; hold on;
    plot(xR(2,:),imag(P_sct),'LineWidth',1.2); 
subplot(1,3,3)
  plot(xR(2,:), abs(Pfsct),'LineWidth',1.2); 
  title('Absolute value of P^{sct}'); axis tight; hold on;
  plot(xR(2,:), abs(P_sct),'LineWidth',1.2); 
  axis tight; hold off
legendtitle1 = sprintf('f=%.2fHz FiniteDiff', f);
legendtitle2 = sprintf('f=%.2fHz Analytic  ', f);
legend(legendtitle1,legendtitle2,'Location','Best');
