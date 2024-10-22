function[outp]=trans(inp);

nx=size(inp,2);
nz=size(inp,1);
outp1(:,1:nx/2) = inp(:,nx/2+1:nx);
outp1(:,nx/2+1:nx) = inp(:,1:nx/2);
outp(1:nz/2,:) = outp1(nz/2+1:nz,:);
outp(nz/2+1:nz,:) = outp1(1:nz/2,:);