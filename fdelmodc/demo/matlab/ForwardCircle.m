function [P_inc,P_sct,xR] = ForwardCircle( xS, xR, gsize, dxR, c_0, c_s, rho_0, rho_s, f )

% xS    = source position
% dxR   = receiver grid size
% c_0   = wave speed in embedding
% c_s   = wave speed of scattering object
% rho_0 = mass density in embedding
% rho_s = mass density of scattering object
% f     = temporal frequency;

wavelength = c_0 / f;                   % wavelength
s          = 1e-16 + 1i*2*pi*f;         % LaPlace parameter
gam_0      = s/c_0;                     % propagation coefficient 0
gam_s      = s/c_s;                     % propagation coefficient s
a          = 40;                        % radius circle cylinder 
depth      = 20;                        % distance betwen circle and origin

disp(['wavelength = ' num2str(wavelength)]);

% Thorbecke coordinates 
   xSource1   = xS(1);                   % source position
   xSource2   = xS(2);
   xReceiver1 = xR(1);                     % receiver positions
   nr=(gsize(2)-1)/2;
   xReceiver2 = (-nr:nr) * dxR;        
   NR         = size(xReceiver2,2); 
   
% Van den Berg coordinates 
   xS(1) = xSource1;                          
   xS(2) = xSource2; 
   
   xR(1,1:NR) = xReceiver1;
   xR(2,1:NR) = xReceiver2(1:NR);

% Compute 2D incident field ---------------------------------------------              
  DIS   = sqrt( (xR(1,:)-xS(1)).^2 + (xR(2,:)-xS(2)).^2 );    
  p_inc =  1/(2*pi).* besselk(0,gam_0*DIS);
 
  % multiply by  rho_0 and plot
  P_inc = rho_0 * p_inc ;                % Take care of factor \rho

 % Make grid
   N1 = gsize(1);                              % number of samples in x_1  
   N2 = gsize(2);                              % number of samples in x_2 
   dx = dxR;                                % with meshsize dx
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
  
   x1 = X1(:,1);   x2 = X2(1,:);
    set(figure(8),'Units','centimeters','Position',[5 5 18 12]); 
    subplot(1,2,1)
       imagesc(x2,x1,cgrid); 
       title('\fontsize{13} c-grid');
       xlabel('x_1 \rightarrow'); 
       ylabel('\leftarrow x_3'); 
       axis('equal','tight');  
       colorbar('hor'); colormap jet;
    subplot(1,2,2)
       imagesc(x2,x1,rhogrid); 
       title('\fontsize{13} \rho-grid');
       xlabel('x_1 \rightarrow'); 
       ylabel('\leftarrow x_3'); 
       axis('equal','tight');  
       colorbar('hor'); colormap jet;



  
  
%-------------------------------------------------------------------------  
% CHECK EXACT INCIDENT FIELD AND BESSEL FUNCTION REPRESENTATION  for r = a
%------------------------------------------------------------------------- 

% Transform Cartesian coordinates to polar ccordinates   
  rS = sqrt(xS(1)^2+xS(2)^2);     phiS = atan2(xS(2),xS(1)); 
  r  = a;                         phi  = 0:.01:2*pi; 
  
% (1) Compute incident wave in closed form --------------------------------
   DIS         = sqrt(rS^2 + r.^2 - 2*rS*r.*cos(phiS-phi));
   p_inc_exact =           1/(2*pi) .* besselk(0,gam_0*DIS);
  dp_inc_exact = - gam_0 *(r-rS*cos(phiS-phi))./DIS ...
                        .* 1/(2*pi) .* besselk(1,gam_0*DIS);

% (2) Compute incident wave as Bessel series with M+1terms --------------
   M = 100;                                 % increase M for more accuracy
  
   p_inc =         besselk(0,gam_0*rS) .* besseli(0,gam_0*r);       
  dp_inc = gam_0 * besselk(0,gam_0*rS) .* besseli(1,gam_0*r); 
  for m = 1 : M;  
       Ib0 = besseli(m,gam_0*r); 
      dIb0 =  gam_0 * (besseli(m+1,gam_0*r) + m/(gam_0*r) * Ib0);
     p_inc =  p_inc + 2 * besselk(m,gam_0*rS)  *  Ib0 .* cos(m*(phiS-phi));                   
    dp_inc = dp_inc + 2 * besselk(m,gam_0*rS) .* dIb0 .* cos(m*(phiS-phi));
  end % m_loop
   p_inc = 1/(2*pi) *  p_inc;
  dp_inc = 1/(2*pi) * dp_inc;
  
% (3) Determine mean error and plot error in domain -----------------------   
  Error_p = p_inc - p_inc_exact;
  disp(['normalized norm of error = '  ...
                  num2str(norm(Error_p(:),1)/norm(p_inc_exact(:),1))]);
  Error_dp = dp_inc - dp_inc_exact;
  disp(['normalized norm of error = '  ...
                  num2str(norm(Error_dp(:),1)/norm(dp_inc_exact(:),1))]);             
                  
  set(figure(9),'Units','centimeters','Position',[5 5 18 14]); 
  subplot(2,1,1)
    angle = phi * 180 / pi;
    semilogy(angle,abs(Error_p)./abs(p_inc_exact));     axis tight;
    xlabel('observation angle in degrees \rightarrow');
    ylabel('abs(p_{inc}-p_{inc}^{exact}) / abs(p_{inc}^{exact}) \rightarrow'); 
  subplot(2,1,2)
    semilogy(angle,abs(Error_dp)./abs(dp_inc_exact));    axis tight;
    title('\fontsize{12} relative error on circle boundary');   
    xlabel('observation angle in degrees \rightarrow');
    ylabel('abs(dp_{inc}-dp_{inc}^{exact}) / abs(dp_{inc}^{exact}) \rightarrow'); 
     
%--------------------------------------------------------------------------  
% COMPUTE SCATTERED FIELD WITH BESSEL FUNCTION REPRESENTATION for r > a
%-------------------------------------------------------------------------- 

% (4) Compute coefficients of series expansion ----------------------------
  Z_0  = c_0 * rho_0;   Z_s  = c_s * rho_s;
  arg0 = gam_0 * a;     args = gam_s *a; 
  A = zeros(1,M+1); 
  for m = 0 : M;                      
    Ib0 = besseli(m,arg0);     dIb0 =  besseli(m+1,arg0) + m/arg0 * Ib0; 
    Ibs = besseli(m,args);     dIbs =  besseli(m+1,args) + m/args * Ibs; 
    Kb0 = besselk(m,arg0);     dKb0 = -besselk(m+1,arg0) + m/arg0 * Kb0;
    A(m+1) = - ((1/Z_s) * dIbs*Ib0 - (1/Z_0) * dIb0*Ibs) ...
              /((1/Z_s) * dIbs*Kb0 - (1/Z_0) * dKb0*Ibs);
  end
  
% (5) Compute scattered field at receivers (data) -------------------------
  rR = sqrt(xR(1,:).^2 + xR(2,:).^2);   phiR = atan2(xR(2,:),xR(1,:));  
  rS = sqrt(xS(1)^2 + xS(2)^2);          phiS = atan2(xS(2),xS(1));
  p_sct = A(1) * besselk(0,gam_0*rS).* besselk(0,gam_0*rR);
  for m = 1 : M;
    factor = 2 * besselk(m,gam_0*rS) .* cos(m*(phiS-phiR));
    p_sct = p_sct + A(m+1) * factor .* besselk(m,gam_0*rR);
  end % m_loop
  p_sct = 1/(2*pi) * p_sct;    
  P_sct = rho_0    * p_sct;                  % Take care of factor \rho
end