function [ P, Vz] = FD_mod( xS, xR, tmod, dtrcv, dx, cgrid, rhogrid, orig)
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
origz=orig(1);
origx=orig(2);
zsrc=xS(2);
xsrc=xS(1);

%write receiver arrays to file
dlmwrite('rcv.txt',xR, ' ');

%compute dt for modeling dt < 0.606*dx/Cmax
Cmax=max(cgrid(:));
dxmod=dx;
dtmax=0.606*dxmod/Cmax;
Sdt=ceil(dtrcv/dtmax);
dtmod=dtrcv/(Sdt);
ntwave=16384;
fmax=0.8/(2*dtrcv); % fmax is 80% of Nyquist frequency
frig=0.75/(2*dtrcv);
flef=0.05/(2*dtrcv);
%tmod=(ntrcv-Sdt+1)*dtrcv;

fileID = fopen('run.scr','w+');
fprintf(fileID,'#!/bin/bash\n');
fprintf(fileID,'export PATH=$HOME/src/OpenSource/bin:/opt/CWP/bin/:.:$PATH\n');
fprintf(fileID,'which fdelmodc\n');
%fprintf(fileID,'set -x\n');
fprintf(fileID,'dt=%e\n',dtmod);
fprintf(fileID,'suaddhead < mod_ro.bin ntrpr=%d ns=%d | \\\n',dims(2), dims(1));
fprintf(fileID,'sushw key=f1,f2,d1,d2,gx,scalco a=%f,%f,%f,%f,%d,-1000 b=0,0,0,0,%d,0 > mod_ro.su\n',origz, origx, dx, dx, int32(origx*1000), int32(dx*1000));
fprintf(fileID,'suaddhead < mod_cp.bin ntrpr=%d ns=%d | \\\n',dims(2), dims(1));
fprintf(fileID,'sushw key=f1,f2,d1,d2,gx,scalco a=%f,%f,%f,%f,%d,-1000 b=0,0,0,0,%d,0 > mod_cp.su\n',origz, origx, dx, dx, int32(origx*1000), int32(dx*1000));
fprintf(fileID,'makewave w=fw fmin=0 flef=%f frig=%f fmax=%f dt=$dt file_out=wavefw.su nt=%d shift=1 scale=0 scfft=1 verbose=1 >& nep\n', flef, frig, fmax, ntwave);
fprintf(fileID,'t0=`grep shift nep | awk ''{print $6}''`\n');
fprintf(fileID,'echo rec_delay for shift in wavelet: t0=$t0\n');
fprintf(fileID,'tmod=$(echo "scale=4; %f+${t0}" | bc -l)\n',tmod);
fprintf(fileID,'export filecp=mod_cp.su\n');
fprintf(fileID,'export filero=mod_ro.su\n');
fprintf(fileID,'export OMP_NUM_THREADS=4\n');
fprintf(fileID,'fdelmodc \\\n');
fprintf(fileID,'file_cp=$filecp file_den=$filero \\\n');
fprintf(fileID,'ischeme=1 \\\n');
fprintf(fileID,'file_src=wavefw.su verbose=1 \\\n');
fprintf(fileID,'dt=$dt \\\n');
fprintf(fileID,'file_rcv=recv.su \\\n');
fprintf(fileID,'rec_type_vz=1 rec_type_p=1 rec_int_vz=2 \\\n');
fprintf(fileID,'rcv_txt=rcv.txt \\\n');
fprintf(fileID,'dtrcv=%e \\\n', dtrcv);
fprintf(fileID,'xsrc=%f zsrc=%f \\\n', xsrc, zsrc);
fprintf(fileID,'src_type=1 tmod=$tmod rec_delay=$t0 \\\n');
fprintf(fileID,'ntaper=100 \\\n');
fprintf(fileID,'left=2 right=2 bottom=2 top=2 \n\n');
fprintf(fileID,'\n');
fprintf(fileID,'sustrip < recv_rp.su > recv_rp.bin\n');
fprintf(fileID,'sustrip < recv_rvz.su > recv_rvz.bin\n');
fprintf(fileID,'surange < recv_rp.su | grep ns | awk ''{print $2}'' > samples\n');
fprintf(fileID,'surange < recv_rp.su | grep traces | awk ''{print $1}'' > traces\n');
fclose(fileID);
!chmod +x run.scr
system('./run.scr');

path = getenv('PATH');
path = [path ':$HOME/src/OpenSource/bin:/opt/CWP/bin/:.:'];
setenv('PATH', path);

% get number of samples and traces
ns=dlmread('samples');
ntr=dlmread('traces');

% Pressure field  
file='recv_rp.bin';
fid=fopen(file,'r'); 
P=fread(fid,[ns,ntr],'float32');
fclose(fid);

% Particle velocity field  
file='recv_rvz.bin';
fid=fopen(file,'r'); 
Vz=fread(fid,[ns,ntr],'float32');
fclose(fid);

end

