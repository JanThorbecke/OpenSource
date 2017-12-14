clear all; clc; close all; clear workspace

% set number of dimensions
global nDIM; nDIM = 2; % set dimension of space

    % set up spatial grid
    N1fd = 1280; N2fd = 1280; dxfd =2.13313822/5;
    x1fd = -(N1fd+1)*dxfd/2 + (1:N1fd)*dxfd;   
    x2fd = -(N2fd+1)*dxfd/2 + (1:N2fd)*dxfd;
    orig = [-(N1fd+1)*dxfd/2,-(N2fd+1)*dxfd/2];
    [X1fd,X2fd] = ndgrid(x1fd,x2fd);   
    
    % load model
%     filn=sprintf('sos_fidmod.bin'); fid = fopen(filn,'r'); sos_fd = fread(fid,[N1fd*N2fd],'float32'); fclose('all');   
%     filn=sprintf('rho_fidmod.bin'); fid = fopen(filn,'r'); rho_fd = fread(fid,[N1fd*N2fd],'float32'); fclose('all');   
%     sos_fd = reshape(sos_fd,[N1fd,N2fd]);
%     rho_fd = reshape(rho_fd,[N1fd,N2fd]);
    sos_fd = ones(N1fd,N2fd)*1500; % speed of sound
    rho_fd = ones(N1fd,N2fd)*1500; % density
    

    % time parameters
    Ntdf = 1024; dtdf = 10^(-3);
 
    % set up acquisition grid
    r_rec = 200; % radius of circle
    Nr = 250; % number of receivers
    rcvr_phi(1:Nr) = (1:Nr) * 2*pi/Nr; % angles
    xR = zeros(2,Nr);
    xR(1,1:Nr) = r_rec * cos(rcvr_phi); 
    xR(2,1:Nr) = r_rec * sin(rcvr_phi);
    xS = xR(:,1) % choose source at first position
    
    % plot the impedance and acquisition geometry
    figure; imagesc(x1fd,x2fd,sos_fd.*rho_fd);
    hold on; scatter(xR(1,:),xR(2,:),'*b');
    hold on; scatter(xS(1,:),xS(2,:),'*r');
    xlabel('x (m)');
    ylabel('y (m)');
    
    % Call finite difference code
    [P, Vz]=FD_mod( xS.', xR.', 0.5, dtdf, dxfd, sos_fd, rho_fd, orig);

    % make a plot
    figure;imagesc(P);
    xlabel('angle (degrees)');
    ylabel('time (t)');
