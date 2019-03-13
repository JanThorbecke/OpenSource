function [ Ptsct, Pfinc, Pfsct, f_out] = FD_mod_grid( xS, xR, ntrcv, dtrcv, dx, cgrid, rhogrid, f )
%Summary of this function goes here
%   Detailed explanation goes here

% save Velocity and density grid
dims = size(cgrid);
fileID =fopen('mod_cp.bin','w+','l');
fwrite(fileID,cgrid,'float32');
fclose(fileID);
fileID =fopen('mod_ro.bin','w+','l');
fwrite(fileID,rhogrid,'float32');
fclose(fileID);

% Compute sizes for makemod
sizez=(dims(1)-1)*dx;
sizex=(dims(2)-1)*dx;
origz=-(dims(1)-1)*dx/2
origx=-(dims(2)-1)*dx/2
zrcv1=xR(1)-origz; 
zrcv2=xR(1)-origz; 
xrcv1=0.0;
xrcv2=xrcv1+(dims(2)-1)*dx; 
zsrc=xS(1)-origz;
xsrc=xS(2)-origx;
tmod=(ntrcv-1)*dtrcv;

%compute dt for modeling dt < 0.606*dx/Cmax
Cmax=max(cgrid(:));
dxmod=dx;
dtmax=0.606*dxmod/Cmax;
dtmod=dtrcv/(ceil(dtrcv/dtmax));
ntwave=16384;
ntfft=4*ntrcv;

fileID = fopen('run.scr','w+');
fprintf(fileID,'#!/bin/bash\n');
fprintf(fileID,'export PATH=$HOME/src/OpenSource/bin:/opt/CWP/bin/:.:$PATH\n');
fprintf(fileID,'which fdelmodc\n');
fprintf(fileID,'which surange\n');
fprintf(fileID,'set -x\n');

fprintf(fileID,'df=0.25\n');
fprintf(fileID,'dt=%e\n',dtmod);
fprintf(fileID,'suaddhead < mod_ro.bin ntrpr=%d ns=%d | \\\n',dims(2), dims(1));
fprintf(fileID,'sushw key=d1,d2,gx,scalco a=%f,%f,%d,-1000 b=0,0,%d,0 > mod_ro.su\n',dx, dx, int32(xrcv1*1000), int32(dx*1000));
fprintf(fileID,'suaddhead < mod_cp.bin ntrpr=%d ns=%d | \\\n',dims(2), dims(1));
fprintf(fileID,'sushw key=d1,d2,gx,scalco a=%f,%f,%d,-1000 b=0,0,%d,0 > mod_cp.su\n',dx, dx, int32(xrcv1*1000), int32(dx*1000));
fprintf(fileID,'makewave w=fw fmin=0 flef=6 frig=94 fmax=100 dt=$dt file_out=wavefw.su nt=%d t0=0.4 scale=0 scfft=1 verbose=1\n', ntwave);
fprintf(fileID,'makewave w=fw fmin=0 flef=6 frig=94 fmax=100 dt=%f file_out=wavefwdt.su nt=%d t0=0.4 scale=0 scfft=1 verbose=1\n', dtrcv, ntfft);
fprintf(fileID,'export filecp=mod_cp.su\n');
fprintf(fileID,'export filero=mod_ro.su\n');
fprintf(fileID,'export OMP_NUM_THREADS=2\n');
fprintf(fileID,'fdelmodc \\\n');
fprintf(fileID,'file_cp=$filecp file_den=$filero \\\n');
fprintf(fileID,'ischeme=1 \\\n');
fprintf(fileID,'file_src=wavefw.su verbose=1 \\\n');
fprintf(fileID,'file_rcv=recvgrid.su \\\n');
fprintf(fileID,'rec_type_vz=0 rec_type_p=1 rec_int_vz=2 \\\n');
fprintf(fileID,'xrcv1=%d xrcv2=%d zrcv1=%d zrcv2=%d \\\n',xrcv1,xrcv2,zrcv1,zrcv2);
fprintf(fileID,'dxrcv=%d \\\n', dx);
fprintf(fileID,'dtrcv=%e \\\n', dtrcv);
fprintf(fileID,'xsrc=%d zsrc=%d\\\n', xsrc, zsrc);
fprintf(fileID,'src_type=1 tmod=%e \\\n', tmod);
fprintf(fileID,'ntaper=100 \\\n');
fprintf(fileID,'left=2 right=2 bottom=2 top=2\n\n');
fprintf(fileID,'\n');
fprintf(fileID,'makemod file_base=hom.su \\\n');
fprintf(fileID,'cp0=1500 ro0=3000 sizex=%d sizez=%d dx=%.1f dz=%.1f orig=0,0 verbose=1\n', sizex,sizez, dxmod,dxmod);
fprintf(fileID,'export filecp=hom_cp.su\n');
fprintf(fileID,'export filero=hom_ro.su\n');
fprintf(fileID,'fdelmodc \\\n');
fprintf(fileID,'file_cp=$filecp file_den=$filero \\\n');
fprintf(fileID,'ischeme=1 \\\n');
fprintf(fileID,'file_src=wavefw.su verbose=1 \\\n');
fprintf(fileID,'file_rcv=recvhom.su \\\n');
fprintf(fileID,'rec_type_vz=0 rec_type_p=1 \\\n');
fprintf(fileID,'xrcv1=%d xrcv2=%d zrcv1=%d zrcv2=%d \\\n',xrcv1,xrcv2,zrcv1,zrcv2);
fprintf(fileID,'dxrcv=%d \\\n', dx);
fprintf(fileID,'dtrcv=%e \\\n', dtrcv);
fprintf(fileID,'xsrc=%d zsrc=%d\\\n', xsrc, zsrc);
fprintf(fileID,'src_type=1 tmod=%e \\\n', tmod);
fprintf(fileID,'ntaper=100 \\\n');
fprintf(fileID,'left=2 right=2 bottom=2 top=2\n\n');
fprintf(fileID,'suop2 recvgrid_rp.su recvhom_rp.su op=diff > recv_rp.su\n');
fprintf(fileID,'sustrip < recv_rp.su > recv_rp.bin\n');
fprintf(fileID,'sustrip < recvhom_rp.su > recvhom_rp.bin\n');
fprintf(fileID,'sustrip < wavefwdt.su > wavefwdt.bin\n');
fclose(fileID);
!chmod +x run.scr
system('./run.scr');

path = getenv('PATH');
path = [path ':$HOME/src/OpenSource/bin:/opt/CWP/bin/:.:'];
setenv('PATH', path);

% Scattered field  
ns=ntrcv;
ntr=dims(2);
file='recv_rp.bin';
fid=fopen(file,'r'); 
Ptsct=fread(fid,[ns,ntr],'float32');
fclose(fid);

% Direct field  
file='recvhom_rp.bin';
fid=fopen(file,'r'); 
Ptinc=fread(fid,[ns,ntr],'float32');
fclose(fid);

df=double(1.0/(ntfft*dtrcv));
a=round(f/df)+1;
f_out=(a-1)*df
fprintf('selected discrete frequency P %f\n', (a-1)*df)
% f2=a*df %these are the selected discrete frequencies
Pfft=fft(Ptsct,ntfft);
Pfsct=Pfft(a,:);
%
Pfft=fft(Ptinc,ntfft);
Pfinc=Pfft(a,:);

ns=ntfft;
ntr=1;
file='wavefwdt.bin';
fid=fopen(file,'r'); 
Wt=fread(fid,[ns,ntr],'float32');
fclose(fid);

df=double(1.0/(ntfft*dtrcv));
c=round(f/df)+1; % select frequencies as close as the one's in Pf
fprintf('selected discrete frequency W %f\n', (c-1)*df)
Wfft=fft(Wt,ntfft);
Wf=Wfft(c,1);

Pfsct = (Pfsct)./(1.0*Wf); % deconvolve for the wavelet
Pfinc = (Pfinc)./(1.0*Wf); % deconvolve for the wavelet

end

