close all
clear all

disp('For generating Green.mp4, choose 1')
disp('For generating Tud.mp4, choose 2')
disp('For generating Wuv.mp4, choose 3')
l=input('Your choice:  ');

axis tight manual
if l==1 
v=VideoWriter('Green.mp4','Motion JPEG AVI');
end
if l==2 
v=VideoWriter('Tud.mp4','Motion JPEG AVI');
end
if l==3 
v=VideoWriter('Wuv.mp4','Motion JPEG AVI');
end

colormap('default');
open(v);

dx	=	0.0002;
nx	=	2048;
nz	=	nx;
nx	=	5/4*nx;
dkx	=	2*pi/dx/nx;
dkz	=	2*pi/dx/nz;
dt	=	50*10^(-6);
cc 	=	[1500, 2500];
resam	=	200;
if l==3
dt	=	25*10^(-6);
cc 	=	[2000, 2000, 1200, 2500, 1400];
resam	=	100;
end
nt	=	length(cc);
	
c	=	zeros(1,resam*nt);

for n=1:nt
for nn=1:resam
c((n-1)*resam+nn)=cc(n);
end
end

nt	=	nt*resam;
ntt	=	nt;
if l==2
ntt=nt-round(resam/8);
end
dt	=	dt/resam;
beta	=	1000*ones(1,nt);
xp	=	0:1:nx/2-1;
zp	=	0:1:nz/2-1;
xm	=	-nx/2:1:-1;
zm	=	-nz/2:1:-1;
xpm	=	dx*[xm xp];
zpm	=	dx*[zm zp];

% avoid division by zero
xp(1)	=	0.05;	
zp(1)	=	0.05;	

kxpm	=	dkx*[xp xm];
kzpm	=	dkz*[zp zm];
[kx,kz]	=	meshgrid(kxpm,kzpm);

% define spatial wavelet
k0	=	2*pi*100;
if l==2
k0	=	2*pi*200;
end
ks	=	kx.*kx+kz.*kz;
k	=	sqrt(ks);
rickerk	=	1/k0^3*ks.*exp(-1/k0^2*ks);	


source	=	zeros(nz,nx);
source(1,1)=1;

if l==2
% define complex source (TUD)
% Letter T:
for p=1:81
source(1,p)=1;
end
for p=1:101
source(p,41)=1;
end
% Letter U:
for p=1:61;
source(60+round(41*sin(p*pi/61)),161+round(41*cos(p*pi/61)))=1;
source(p,120)=1;
source(p,202)=1;
end
% Letter D:
for p=1:101
source(p,250)=1;
source(51+round(50*cos(p*pi/101)),290+round(50*sin(p*pi/101)))=1;
end
for p=1:41
source(1,250+p)=1;
source(101,250+p)=1;
end
end

sourcek	=	fft2(source);	

% initialize propagator matrix (equation 11.7) and multiply with source
	Wuuk	=	ones(nz,nx).*rickerk.*sourcek;
	Wuvk	=	zeros(nz,nx);
	Wvuk	=	zeros(nz,nx);
	Wvvk	=	ones(nz,nx).*rickerk.*sourcek;
if l==3
	Wuuk	=	cos(-c(1)*resam*dt*k).*rickerk;
	Wuvk	=	-1/(beta(1)*c(1))*sin(-c(1)*resam*dt*k)./k.*rickerk;
	Wvuk	=	beta(1)*c(1)*sin(-c(1)*resam*dt*k).*k.*rickerk;
	Wvvk	=	Wuuk;
end

clims	=	1.0*10^(-16)*[-1 1];


for n=1:ntt
	Wuvx	=	real(trans(ifft2(Wuvk)));	
	Wuv(:,n)=	Wuvx(:,nx/2+1);
if l==1
% extract snapshot of Green's function (equation 11.22)
	imagesc(xpm,zpm,-Wuvx,clims/sqrt(n))
elseif l==2
% extract snapshot of Green's function, convolved with complex source function
	imagesc(xpm,zpm,-Wuvx,clims)
else
	imagesc(xpm,zpm,Wuvx,clims/sqrt(n))
end
	frame = getframe(gcf);
	writeVideo(v,frame);
% define propagator for time-slab n (equation B1)
	wuuk	=	cos(c(n)*dt*k);
	wuvk	=	-1/(beta(n)*c(n))*sin(c(n)*dt*k)./k;
	wvuk	=	beta(n)*c(n)*sin(c(n)*dt*k).*k;
	wvvk	=	wuuk;

	Wuukh	=	Wuuk;
	Wuvkh	=	Wuvk;
	Wvukh	=	Wvuk;
	Wvvkh	=	Wvvk;
% apply one recursion step (equation 11.10)
	Wuuk	=	wuuk.*Wuukh + wuvk.*Wvukh;
	Wuvk	=	wuuk.*Wuvkh + wuvk.*Wvvkh;
	Wvuk	=	wvuk.*Wuukh + wvvk.*Wvukh;
	Wvvk	=	wvuk.*Wuvkh + wvvk.*Wvvkh;
end

close(v);

if l==1
% plot Green's function (equation 11.22) as z,t diagram 
figure(2)
wiggle(-Wuv(:,1:10:nt),[dx,10*dt],[-nz/2*dx,0],5,1,'k-')
print -deps figure2.eps
end

if l==3
% plot component Wuv function as z,t diagram
figure(3)
wiggle(Wuv(:,1:10:nt),[dx,10*dt],[-nz/2*dx,0],5,1,'k-')
print -deps figure3.eps
% extract and plot Green's function (equation 11.22) as z,t diagram 
Wuv(:,1:resam)=zeros(nz,resam);
figure(4)
wiggle(-Wuv(:,1:10:nt),[dx,10*dt],[-nz/2*dx,0],5,1,'k-')
print -deps figure4.eps
end
