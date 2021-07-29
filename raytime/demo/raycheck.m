clear; close all; clc;

nx = 801;
sx = 4000;
nt = 1634;
dx = 25;
dt = 0.003;
xvec = sx:dx:sx+dx*(nx-1);
tvec = 0:dt:(nt-1)*dt;

figure;
FD = fopen('FD_rp.bin');
gm = reshape(fread(FD,'float32'),nt,nx);
fclose(FD);
imagesc(xvec,tvec,gm);
colormap gray
ylim([0 5])
caxis([-10 10])
set(gca,'FontSize',18)
set(gca,'XAxisLocation','top');
xlabel('x [m]');
ylabel('t [s]')
title('Saga depth 4500')

hold on;
FD = fopen('pray0_time.bin');
g2 = reshape(fread(FD,'float32'),nx,1);
fclose(FD);
plot(xvec,g2,'--k','LineWidth',2);

hold on;
FD = fopen('Jray0_time.bin');
g2 = reshape(fread(FD,'float32'),nx,1);
fclose(FD);
plot(xvec,g2,'--g','LineWidth',2);

hold on;
FD = fopen('Jray0_timeStandard.bin');
g2 = reshape(fread(FD,'float32'),nx,1);
fclose(FD);
plot(xvec,g2,'--c','LineWidth',2);

legend('FD method','Jesper with 0.5*T2', 'Standard Jesper')

% 
% nx = 801;
% sx = 0;
% nt = 2049;
% dx = 25;
% dt = 0.004;
% xvec = sx:dx:sx+dx*(nx-1);
% tvec = 0:dt:(nt-1)*dt;
% 
% figure;
% FD = fopen('saga_shot.bin');
% gm = reshape(fread(FD,'float32'),nt,nx);
% fclose(FD);
% imagesc(xvec,tvec,gm);
% colormap french
% ylim([1.7 4])
% caxis([-0.00001 0.00001])
% set(gca,'FontSize',18)
% set(gca,'XAxisLocation','top');
% xlabel('x [m]');
% ylabel('t [s]')
% title('SAGA')
% 
% hold on;
% FD = fopen('saga_fd.bin');
% g2 = reshape(fread(FD,'float32'),nx,1);
% fclose(FD);
% plot(xvec,g2,'--k','LineWidth',2);
% 
% hold on;
% FD = fopen('saga_T1.bin');
% g2 = reshape(fread(FD,'float32'),nx,1);
% fclose(FD);
% plot(xvec,g2,'--g','LineWidth',2);
% 
% hold on;
% FD = fopen('saga_T2.bin');
% g2 = reshape(fread(FD,'float32'),nx,1);
% fclose(FD);
% plot(xvec,g2,'--c','LineWidth',2);
% 
% legend('old','T1','T2')