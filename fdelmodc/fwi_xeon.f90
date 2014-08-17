! Module containing different subroutines for Time Domain Full Waveform Inversion.
! Also include 2-4 stagerred grid finite difference routines of both P and SH modelling
! where the pml absorbing conditions are modified from Chaiwoot's CSIM code
! Gradient calculations can be done for different obj functions and also different inversion schemes can
! be used..
! Author: Pawan Bharadwaj
!	  p.b.pisupati@tudelft.nl
!         Department of Geotechnology
!         Applied Geophysics and Petrophysics
!         CiTG - Delft University of Technology
!         
! Date  : October, 2012

module fwi_xeon
use data_global
use my_func
public
integer,public		:: nx,nz,iband
! nz and nx	--	total size of the model
! iband		-- 	# freq band which is being inverted
integer,private 	:: na_pml, na_pad, &
  na_new, na_store,nx0,nz0
real, private 		:: dtx
! dtx		=	dt*dx
integer, private 	:: ntaper
integer,PRIVATE 	:: nxm1,nxm2,nzm1,nzm2,nxm10,nxm20,nzm10,nzm20 !model dimensions
integer, private 	::nzm,nxm
parameter(ntaper=10) !need not change this one .. 
parameter(na_pml=40) !grids added for pmls
parameter(na_pad=2)  !grids added for padding
parameter(na_new=na_pml+na_pad) 
! na_new 	--	total number of grids added to the model of PML+padding
real, dimension(:,:), allocatable, PRIVATE &
			:: c1,c2,cl,den, &
			   cl41,cl42, &
			   exL,exR, &
			   eyT,eyB 
! ex* and ey*	--	damping constants for the pml conditions 
!$OMP THREADPRIVATE(c1,c2,cl,den,cl41,cl42,exL,exR,eyT,eyB)
! tunneling parameters
integer, private		:: ntundep,ntunrad,ntungap
! c41 and c42	--	constants for stagerred 2D stagerred grid kernel
real, private		:: c41, c42
parameter(c41=9./8.)
parameter(c42=-1./24.)
!ILLUMIATIONS
real, private,allocatable,dimension(:,:) &
			:: illu_src,illu_rec
real, private,allocatable,dimension(:,:,:) &
			:: illu

contains

subroutine do_fw_modelling(flag)
! This subroutine does forward modelling for istage
!
! Author : Pawan Bharadwaj
!	   p.b.pisupati@tudelft.nl
!
! Input:
!	flag	-1  -- uses original velocity and density model for modelling/
!    	  	 0  -- uses c (velocity) and rho (density) for modelling

! Output: global variables
!	csgs_cal    -- common shot gathers if flag.eq.0
! 	csgs_obs    -- common shot gathers if flag.eq.-1

use data_global
use omp_lib
implicit none
		integer 			:: is,istg
		integer,intent(in) 		:: flag
		character(len=200) 		:: junk
		real, dimension(nt,ng,ns)	:: csgs_temp	
		real, dimension(nt)		:: src_temp
csgs_temp(:,:,:)=0.0
src_temp(:)=0.0
if(my_rank.eq.0) then
!	write(*,*)'forward modelling ...'
end if
do istg=1,nstage/nump
	istage=1+my_rank+(istg-1)*nstage/nump
	!$OMP PARALLEL NUM_THREADS(num_threads) DEFAULT(SHARED) &
	!$OMP PRIVATE(junk,src_temp)
	!$OMP DO
	do is=ns/nstage*(istage-1)+1,ns/nstage*istage
		if(flag.eq.-1) then
			write(*,*)'source position(x,z)=',xs(is),'m,',zs(is),'m'
		end if
	        if(iband.ne.0.and.flag.eq.0) then
        		call do_bw_filter1d(src_temp,src_sign,bands_freq(iband,1),bands_freq(iband,2),dt)
       	 	elseif(flag.eq.-1) then
			src_temp=src_sign
		endif
		call forw(vel_flag=flag,csg=csgs_temp(:,:,is),is=is,source_wav=src_temp)
	if(flag.eq.-1) then
		call filename(junk,trim(out_dir)//'/csg',is,'')
		call makesu2d(csgs_temp(:,:,is),junk,0.0,0.0,dx,dx,ng,nt,OMP_get_thread_num()+10)
	endif
	enddo
	!$OMP END DO
	!$OMP END PARALLEL
enddo
call mpi_add_3d(csgs_temp)
if(flag.eq.-1) then
	csgs_obs=csgs_temp;
else if(flag.eq.0) then
	csgs_cal=csgs_temp
endif
end subroutine do_fw_modelling

subroutine vel_den_shot(stage,flag)
use data_global
use omp_lib
implicit none
!x coordinate position of the source -- x_pos
!flag= -1  -- uses original velocity and density model for modelling/
!    =  0  -- uses c (velocity) and rho (density) for modelling
integer, intent(in) 	:: stage,flag
integer 		:: ix,iz,izz,ixx
real 			:: factor_att, const,tmp,alpha,alpha_max0,tmp0
character(len=200) 	:: junk
dtx=dt/dx
na_store=na_pml
nx0=nx_in+2*na_new
nz0=nz_in+2*na_new
allocate	(c1(nx0,nz0),c2(nx0,nz0),&
		cl(nx0,nz0),den(nx0,nz0), &
		cl41(nx0,nz0),cl42(nx0,nz0), &
		exL(na_pml,nz0),exR(na_pml,nz0), &
		eyT(nx0,na_pml),eyB(nx0,na_pml))
nxm=nx_in
nzm=nz_in
nxm1=na_new+1
nzm1=na_new+1
nxm2=nxm1+nxm-1
nzm2=nzm1+nzm-1
nxm10=nxm1-na_pad
nzm10=nzm1-na_pad
nxm20=nxm2+na_pad
nzm20=nzm2+na_pad
nx=nxm+2*na_new
nz=nzm+2*na_new
c1(:,:)=0.0
den(:,:)=0.0
if(flag.eq.-1) then
	c1(na_new+1:nx_in-na_new,na_new+1:nz_in-na_new)=transpose(velocity)
	den(na_new+1:nx_in-na_new,na_new+1:nz_in-na_new)=transpose(density)
else if(flag.eq.0) then
	c1(na_new+1:nx_in-na_new,na_new+1:nz_in-na_new)=transpose(reshape&
                        (c(nxm*nzm*(stage-1)+1:nxm*nzm*(stage)),(/nz_in,nx_in/)))
	den(na_new+1:nx_in-na_new,na_new+1:nz_in-na_new)=transpose(reshape&
                        (rho(nxm*nzm*(stage-1)+1:nxm*nzm*(stage)),(/nz_in,nx_in/)))
end if
!padding velocity and density layers
do ix=1,nx0
	do iz=1,na_new
		c1(ix,iz)=c1(ix,na_new+1)
		c1(ix,nz0-iz+1)=c1(ix,na_new+nz_in)
		den(ix,iz)=den(ix,na_new+1)
		den(ix,nz0-iz+1)=den(ix,na_new+nz_in)
	end do
end do
do iz=1,nz0
	do ix=1,na_new
		c1(ix,iz)=c1(na_new+1,iz)
		c1(nx_in+na_new+ix,iz)=c1(na_new+nx_in,iz)
		den(ix,iz)=den(na_new+1,iz)
		den(nx_in+na_new+ix,iz)=den(na_new+nx_in,iz)
	end do
end do

if(free.eq.1) then
	c1(:,1:na_new)=0;
	den(:,1:na_new)=200; !density of air in kg/m3
end if
!inserting the tunnel in actual density and velocity models
	ntundep=int(tundep/dx)
	ntunrad=int(tunrad/dx)
	ntungap=int(tungap/dx)
if(tundep.ne.0.and.flag.eq.-1) then
	do ix=1,na_new+nxs(1+ns/nstage*(stage-1))-ntungap
		do iz=na_new+ntundep,na_new+ntundep+ntunrad
			c1(ix,iz)=tunvel
			den(ix,iz)=tunrho
		end do
	end do
end if 

!$OMP MASTER
call makesu2d(transpose(c1)&
          ,trim(out_dir)//'/velocity_used',0.0,0.0,dx,dx,nx0,nz0,OMP_get_thread_num()+10)
!$OMP END MASTER
! ! IMPORTANT **** CONSTANTS FOR P or SH case
if(wtd.eq.1) then
! need to change density - but velocity will be taken as
! S-wave velocity
	do ix=1,nx0
		do iz=1,nz0
			if(c1(ix,iz).ne.0.and.den(ix,iz).ne.0) then
				den(ix,iz)=1/c1(ix,iz)/c1(ix,iz)/den(ix,iz) !fake density
			end if
		end do
	end do
end if
do ix=1,nx
	do iz=1,nz
		cl(ix,iz)=dtx/den(ix,iz)
		cl41(ix,iz)=cl(ix,iz)*c41
		cl42(ix,iz)=cl(ix,iz)*c42
	enddo
enddo
do iz=1,nz
	do ix=1,nx
		c2(ix,iz)= (c1(ix,iz)**2)*den(ix,iz)*dtx
	enddo
enddo
!   exL
factor_att=0.1
const=2.302585*dtx
alpha_max0=sqrt(3.0)*vmax*(8.0/15.0-0.03*na_pml &
             +na_pml*na_pml/1500.0)*dtx
do ix=1,na_pml
	tmp=(na_pml+1-ix)/float(na_pml+1)
	tmp=tmp*tmp
do iz=1,nz
	alpha=alpha_max0*tmp
	exL(ix,iz)=1.0-alpha
enddo
enddo
!  exR
do ix=nxm20+1,nx
	ixx=nx-ix+1
	tmp=(na_pml+1-ixx)/float(na_pml+1)
	tmp=tmp*tmp
	do iz=1,nz
		alpha=alpha_max0*tmp
		exR(ixx,iz)=1.0-alpha
	enddo
enddo
!  eyT
do iz=1,na_pml
	tmp0=(na_pml+1.0-iz)/float(na_pml+1)
	tmp0=tmp0*tmp0
	do ix=1,nx
		alpha=alpha_max0*tmp0
		eyT(ix,iz)=1.0-alpha
	enddo
enddo
! ! c  eyB
do iz=nzm20+1,nz
	izz=nz-iz+1
	tmp0=(na_pml+1.0-izz)/float(na_pml+1)
	tmp0=tmp0*tmp0
	do ix=1,nx
		alpha=alpha_max0*tmp0
		eyB(ix,izz)=1.0-alpha
	enddo
enddo
end subroutine vel_den_shot

subroutine forw(vel_flag,csg,is,ig,source_wav,bord_UD,bord_RL,snaps,illu)
!vel_flag = -1 uses velocity and density
!vel_flag =  0 uses c and rho for modelling
! if 'is' is present .. source at nxs(is) and nzs(is)
! if 'ig' and is are present then source at nxg(ig,is) and nzg(ig,is)
! if is <0 .. all sources blasted simultaneoulsy - for calculating illumination
! if ig <0 .. all receivers of is are blasted simultaneoulsy - for calculating illumination
use data_global
use omp_lib
! x_pos and z_pos are coordinate postions of source in the original model
implicit none
integer, intent(in)		:: vel_flag
integer :: it,ix,i
real, dimension(nt,nx_in,nz_in),intent(out)&
				::snaps
real, dimension(nz_in,nx_in),intent(inout) &
				::illu
integer,intent(in)		:: is,ig ! index of the source in nxs variable
real, dimension(nt),intent(in)	::source_wav
real, dimension(:,:), allocatable&
				:: p2xL,p2xR,p2yL,p2yR,&
                            	   p2xT,p2yT,p2xB,p2yB,p2,u1,w1,p2prev
real, dimension(nt,nz_in,4+nx_in), intent(out) &
				:: bord_UD
real, dimension(nt,nx_in,4+nz_in), intent(out) &
				:: bord_RL
real, dimension(nt,ng), intent(out)&
				:: csg
optional 			:: bord_UD,bord_RL,snaps,csg,ig,is,illu
integer k
real, dimension(int(ns/nstage)) :: src_array_taper
real, dimension(ng) :: rec_array_taper

call taper(ns/nstage,sample=src_array_taper,percent=10.0);
call taper(ng,sample=rec_array_taper,percent=10.0);
nx0=nx_in+2*na_new
nz0=nz_in+2*na_new
nx=nx0
nz=nz0
allocate(p2(1:nx0,1:nz0),u1(1:nx0,1:nz0),&
          w1(1:nx0,1:nz0),p2xL(1:na_pml,1:nz0),&
          p2xR(1:na_pml,1:nz0),p2yL(1:na_pml,1:nz0),&
          p2yR(1:na_pml,1:nz0),p2xT(1:nx0,1:na_pml),&
          p2yT(1:nx0,1:na_pml),&
          p2xB(1:nx0,1:na_pml),p2yB(1:nx0,1:na_pml))
allocate(p2prev(1:nx0,1:nz0))
call vel_den_shot(int((is-1)/(ns/nstage))+1,vel_flag)

!     ----------------------------------------------------------------
!    Forward modeling to get synthetic seismograms
!    -----------------------------------------------------------------
 p2(1:nx0,1:nz0)=0.0; p2prev=p2
 u1(1:nx0,1:nz0)=0.0
 w1(1:nx0,1:nz0)=0.0
 p2xL(1:na_pml,1:nz0)=0.0
 p2yL(1:na_pml,1:nz0)=0.0
 p2xR(1:na_pml,1:nz0)=0.0
 p2yR(1:na_pml,1:nz0)=0.0
 p2xT(1:nx0,1:na_pml)=0.0
 p2yT(1:nx0,1:na_pml)=0.0
 p2xB(1:nx0,1:na_pml)=0.0
 p2yB(1:nx0,1:na_pml)=0.0

!  ! !  Main Loop
 dt_loop:do it=1,nt
 if(maxval(p2).gt.10**6) then
 write(*,*)'CHECK STABILITY -- VALUES EXPLODING > 10**6'
 STOP
 end if
!!! adding source at source or receiver position 
if(present(ig).eqv..false.) then
        p2(na_new+nxs(is),na_new+nzs(is))=p2(na_new+nxs(is),na_new+nzs(is))+&
        source_wav(it)!*src_array_taper(is-ns/nstage*(istage-1))
else
         p2(na_new+nxg(ig,is),na_new+nzg(ig,is))=p2(na_new+nxg(ig,is),na_new+nzg(ig,is))+&
        source_wav(it)
end if
 !boundary conditions for the tunnel wall and free surface (SH WAVE CASE)
 !   u1 --> tauyx
 !   w1 --> tauyz
 !   p  --> vy
 if(tundep.ne.0 .and. wtd.eq.1) then
   w1(1:int(na_new+nxs(is)-ntungap)+1,int(na_new+ntundep-1))=0
   w1(1:int(na_new+nxs(is)-ntungap)+1,int(na_new+ntundep+ntunrad))=0
   u1(int(na_new+nxs(is)-ntungap),&
 int(na_new+ntundep-1):int(na_new+ntundep+ntunrad+1))=0
 end if
 if(free.eq.1 .and. wtd.eq.1) then
 w1(1:nx0,na_new+1)=0
 end if
 !boundary conditions for the free surface (P WAVE CASE)
 if(free.eq.1 .and. wtd.eq.2)then
          do ix=1,nx
            p2(ix,na_new+1)=0.0
          enddo
 endif
 
p2prev=p2
 call move(p2,u1,w1, &
          p2xL, &
          p2xR, &
          p2yL, &
          p2yR, &
          p2xT, &
          p2yT, &
          p2xB, &
          p2yB)
!$OMP MASTER
if(it.eq.int(nt/3)) then
call makesu2d(transpose(p2(na_new+1:na_new+nx_in,na_new+1:na_new+nz_in)),&
             trim(out_dir)//'/snap_forw',0.0,0.0,dx,dx,nx_in,nz_in,OMP_get_thread_num()+10)
end if
!$OMP END MASTER
 if(present(snaps))then
 
 snaps(it,:,:)=p2(na_new+1:nx_in+na_new,na_new+1:na_new+nz_in)
 
 end if
if(present(illu)) then
illu=illu+(transpose((p2(na_new+1:nx_in+na_new,na_new+1:na_new+nz_in)&
          -p2prev(na_new+1:nx_in+na_new,na_new+1:na_new+nz_in)/dt)**2))
end if
! ! 
! ! !c
! !         if(nfree.eq.1)then
! !           do ix=1,nx
! !              iz=nsurf(ix)-1
! !              w1(ix,iz) = w1(ix,nsurf(ix)+1)
! ! ! c              w1(ix,iz)=w1(ix,iz)+
! ! ! c     1           cl41(ix,iz)*(-p2(ix,iz))
! ! ! c     2          +cl42(ix,iz)*(p2(ix,iz+2)+p2(ix,iz+3))
! !            enddo
! !         end if
! !c
! ! 	if(ispre.ne.1)then
! !
! !  output the seismogram (pressure)
! !
 if(present(csg)) then
 do k=1,ng
 csg(it,k)=p2(na_new+nxg(k,is),na_new+nzg(k,is))!*rec_array_taper(k)
 enddo
 end if
! ! 	endif
! 
       if(present(bord_UD).and.present(bord_RL)) then
 !
 !      store records and boundary value
 bord_UD(it,:,1)=p2(na_new+1        ,na_new+1:na_new+nz_in)
 bord_UD(it,:,4)=p2(na_new+nx_in    ,na_new+1:na_new+nz_in)
 bord_UD(it,:,2)=p2(na_new+2        ,na_new+1:na_new+nz_in)
 bord_UD(it,:,3)=p2(na_new+nx_in-1  ,na_new+1:na_new+nz_in)
 if(it.eq.nt-1) then
 do k=1,nx_in
 bord_UD(1,:,4+k)=p2(na_new+k,na_new+1:na_new+nz_in)
 end do
 end if
 bord_RL(it,:,1)=p2(na_new+1:na_new+nx_in  ,      na_new+1)
 bord_RL(it,:,4)=p2(na_new+1:na_new+nx_in  ,  na_new+nz_in)
 bord_RL(it,:,2)=p2(na_new+1:na_new+nx_in  ,      na_new+2)
 bord_RL(it,:,3)=p2(na_new+1:na_new+nx_in  ,na_new+nz_in-1)
 if(it.eq.nt) then
 do k=1,nz_in
 bord_RL(1,:,4+k)=p2(na_new+1:na_new+nx_in,na_new+k)
 end do
 end if
         endif
 !  End of Main Loop
 end do dt_loop

deallocate(c1,c2,cl,den, &
       cl41,cl42,&
       exL,exR, &
       eyT,eyB)


end subroutine forw




subroutine back(vel_flag,grad,snaps,is1,csg1,is2,csg2,illu)
!vel_flag = -1 uses velocity and density
!vel_flag =  0 uses c and rho for modelling
use data_global
use omp_lib
! x_pos and z_pos are coordinate postions of source in the original model
implicit none
integer :: ix,iz,it,k
integer, intent(in) :: vel_flag
real, dimension(nz_in,nx_in),intent(inout) :: grad
real, dimension(nz_in,nx_in) :: gradtemp
integer :: x_pos1,z_pos1,x_pos2,z_pos2
real, dimension(nt,nx_in,nz_in),intent(in)::snaps
real, dimension(nt,nx_in,nz_in)              ::dsnaps
integer,intent(in)::is1,is2 !index number of the sources
real, dimension(nt,ng),intent(in):: csg1,csg2
real, dimension(nt,ng)           :: csg
real, dimension(nz_in,nx_in),intent(inout) &
				::illu
! is2 is zero then migrate shot record at is1
! if is2 is not equal to zero then migrate shotrecord(is1-is2)
real, dimension(:,:), allocatable :: p2xL,p2xR,p2yL,p2yR,&
                            p2xT,p2yT,p2xB,p2yB,p2,u2,w2 ! for back propagating boundary data
real, dimension(:,:), allocatable :: p1xL,p1xR,p1yL,p1yR,&
                            p1xT,p1yT,p1xB,p1yB,p1,u1,w1,p1prev ! for back propagating recorded data
optional :: is2,csg2,snaps,illu,grad

nx0=nx_in+2*na_new
nz0=nz_in+2*na_new
nx=nx0
nz=nz0

if(present(is2).and.present(csg2)) then
csg=csg1-csg2
x_pos2=nxs(is2)
z_pos2=nzs(is2)
else
csg=csg1
end if

x_pos1=nxs(is1)
z_pos1=nzs(is1)

allocate(p1(1:nx0,1:nz0),u1(1:nx0,1:nz0),&
          w1(1:nx0,1:nz0),p1xL(1:na_pml,1:nz0),&
          p1xR(1:na_pml,1:nz0),p1yL(1:na_pml,1:nz0),&
          p1yR(1:na_pml,1:nz0),p1xT(1:nx0,1:na_pml),&
          p1yT(1:nx0,1:na_pml),&
          p1xB(1:nx0,1:na_pml),p1yB(1:nx0,1:na_pml))
allocate(p1prev(1:nx0,1:nz0))
call vel_den_shot(int((is1-1)/(ns/nstage))+1,vel_flag) ! getting velocity and density models

!     ----------------------------------------------------------------
!  Initializing values to zero
!    -----------------------------------------------------------------
grad(:,:)=0.0
p1(1:nx0,1:nz0)=0.0; p1prev=p1
u1(1:nx0,1:nz0)=0.0
w1(1:nx0,1:nz0)=0.0
p1xL(1:na_pml,1:nz0)=0.0
p1yL(1:na_pml,1:nz0)=0.0
p1xR(1:na_pml,1:nz0)=0.0
p1yR(1:na_pml,1:nz0)=0.0
p1xT(1:nx0,1:na_pml)=0.0
p1yT(1:nx0,1:na_pml)=0.0
p1xB(1:nx0,1:na_pml)=0.0
p1yB(1:nx0,1:na_pml)=0.0

! time differential of the forward propagated wavefield . . .
if(present(snaps)) then
	do it=1,nt-1
		dsnaps(it,:,:)=(snaps(it+1,:,:)-snaps(it,:,:))/dt
	end do
end if

dt_loop:do it=nt,1,-1 ! backward marching

! input data at the receivers
do k=1,ng
	p1(na_new+nxg(k,is1),na_new+nzg(k,is1))= &
		csg(it,k)+p1(na_new+nxg(k,is1),na_new+nzg(k,is1))
end do

if(present(snaps)) then
	forall(ix=1:nx_in,iz=1:nz_in)
!go_back
! the formula derived from Jerry's book after changing slowness derivative to 
! velocity derivateive.. i.e., glsc from glss
! glsc		--	Gradient of Leastsquare functional wrt 'c'
! glsc		--	Gradient of    ''             ''   wrt 's'
	gradtemp(iz,ix)=-2.0*&!/(c1(na_new+ix,na_new+iz)**3)*&!*den(na_new+ix,na_new+iz))* &
! derivative in time 
       		((p1(na_new+ix,na_new+iz)-p1prev(na_new+ix,na_new+iz))/dt)&
		*dsnaps(it,ix,iz) !gradient result
	end forall
if(tundep.ne.0) then
gradtemp(:,1:int(x_pos1-ntungap)+1)=0;
endif

grad=grad+gradtemp
end if
!------------------------------------------------------------------
! Back propagating data/ data residuals
!------------------------------------------------------------------
if(maxval(p1).gt.10**6) then
	write(*,*)'CHECK STABILITY -- VALUES EXPLODING > 10**6'
	STOP
end if

!boundary conditions for the tunnel wall and free surface (SH WAVE CASE)
!   u1 --> tauyx
!   w1 --> tauyz
!   p  --> vy
if(tundep.ne.0 .and. wtd.eq.1) then
	w1(1:int(na_new+x_pos1-ntungap)+1,int(na_new+ntundep-1))=0
	w1(1:int(na_new+x_pos1-ntungap+1),int(na_new+ntundep+ntunrad))=0
	u1(int(na_new+x_pos1-ntungap),&
		int(na_new+ntundep-1):int(na_new+ntundep+ntunrad+1))=0
end if
if(free.eq.1 .and. wtd.eq.1) then
	w1(1:nx0,na_new+1)=0
end if


!boundary conditions for the free surface (P WAVE CASE)
if(free.eq.1 .and. wtd.eq.2)then
	do ix=1,nx
		p1(ix,na_new+1)=0.0
	enddo
endif

p1prev=p1
call move(p1,u1,w1, &
         p1xL, &
         p1xR, &
         p1yL, &
         p1yR, &
         p1xT, &
         p1yT, &
         p1xB, &
         p1yB)

if(present(illu)) then
	illu=illu+(transpose((p1(na_new+1:nx_in+na_new,na_new+1:na_new+nz_in)&
		-p1prev(na_new+1:nx_in+na_new,na_new+1:na_new+nz_in)/dt)**2))
endif

!$OMP MASTER
if(it.eq.int(nt/3)) then
call makesu2d(transpose(p1(na_new+1:na_new+nx_in,na_new+1:na_new+nz_in)),&
              trim(out_dir)//'/snap_back',0.0,0.0,dx,dx,nx_in,nz_in,OMP_get_thread_num()+10)
end if
!$OMP END MASTER

end do dt_loop

!  End of Main Loop
!=======================================================
deallocate(c1,c2,cl,den, &
       cl41,cl42,&
       exL,exR, &
       eyT,eyB)

end subroutine back

subroutine back_testing(vel_flag,image,bordUD,bordRL,is,csg)
!vel_flag = -1 uses velocity and density
!vel_flag =  0 uses c and rho for modelling
use data_global
use omp_lib
! x_pos and z_pos are coordinate postions of source in the original model
implicit none
integer :: ix,iz,it,k
integer, intent(in) :: vel_flag
real, dimension(nz_in,nx_in),intent(inout) :: image
real, dimension(nz_in,nx_in) :: imagetemp
integer :: x_pos,z_pos
integer,intent(in)::is !index number of the sources
real, dimension(nt,ng),intent(in):: csg
real, dimension(nt,nz_in,4+nx_in), intent(in) :: bordUD
real, dimension(nt,nx_in,4+nz_in), intent(in) :: bordRL
real, dimension(:,:), allocatable :: p2xL,p2xR,p2yL,p2yR,&
                            p2xT,p2yT,p2xB,p2yB,p2,u2,w2,p2prev ! for back propagating boundary data
real, dimension(:,:), allocatable :: p1xL,p1xR,p1yL,p1yR,&
                            p1xT,p1yT,p1xB,p1yB,p1,u1,w1,p1prev ! for back propagating recorded data
nx0=nx_in+2*na_new
nz0=nz_in+2*na_new
nx=nx0
nz=nz0


x_pos=nxs(is)
z_pos=nzs(is)

allocate(p2(1:nx0,1:nz0),u2(1:nx0,1:nz0),&
          w2(1:nx0,1:nz0),p2xL(1:na_pml,1:nz0),&
          p2xR(1:na_pml,1:nz0),p2yL(1:na_pml,1:nz0),&
          p2yR(1:na_pml,1:nz0),p2xT(1:nx0,1:na_pml),&
          p2yT(1:nx0,1:na_pml),&
          p2xB(1:nx0,1:na_pml),p2yB(1:nx0,1:na_pml))
allocate(p2prev(1:nx0,1:nz0))
allocate(p1(1:nx0,1:nz0),u1(1:nx0,1:nz0),&
          w1(1:nx0,1:nz0),p1xL(1:na_pml,1:nz0),&
          p1xR(1:na_pml,1:nz0),p1yL(1:na_pml,1:nz0),&
          p1yR(1:na_pml,1:nz0),p1xT(1:nx0,1:na_pml),&
          p1yT(1:nx0,1:na_pml),&
          p1xB(1:nx0,1:na_pml),p1yB(1:nx0,1:na_pml))
allocate(p1prev(1:nx0,1:nz0))

call vel_den_shot(int((is-1)/(ns/nstage))+1,vel_flag) ! getting velocity and density models

!     ----------------------------------------------------------------
!  Initializing values to zero
!    -----------------------------------------------------------------
imagetemp(:,:)=0.0
image(:,:)=0.0
p2(1:nx0,1:nz0)=0.0; p2prev=p2
u2(1:nx0,1:nz0)=0.0
w2(1:nx0,1:nz0)=0.0
p2xL(1:na_pml,1:nz0)=0.0
p2yL(1:na_pml,1:nz0)=0.0
p2xR(1:na_pml,1:nz0)=0.0
p2yR(1:na_pml,1:nz0)=0.0
p2xT(1:nx0,1:na_pml)=0.0
p2yT(1:nx0,1:na_pml)=0.0
p2xB(1:nx0,1:na_pml)=0.0
p2yB(1:nx0,1:na_pml)=0.0
p1(1:nx0,1:nz0)=0.0; p1prev=p1
u1(1:nx0,1:nz0)=0.0
w1(1:nx0,1:nz0)=0.0
p1xL(1:na_pml,1:nz0)=0.0
p1yL(1:na_pml,1:nz0)=0.0
p1xR(1:na_pml,1:nz0)=0.0
p1yR(1:na_pml,1:nz0)=0.0
p1xT(1:nx0,1:na_pml)=0.0
p1yT(1:nx0,1:na_pml)=0.0
p1xB(1:nx0,1:na_pml)=0.0
p1yB(1:nx0,1:na_pml)=0.0
!write(*,*)'Calculating gradient for shot# =',is1

dt_loop:do it=nt,1,-1 ! backward marching
forall(ix=1:nx_in,iz=1:nz_in)
imagetemp(iz,ix)=-2.0*&!/(c1(na_new+ix,na_new+iz)**3)*&!*den(na_new+ix,na_new+iz))* &
! derivative in time 
       ((p1(na_new+ix,na_new+iz)-p1prev(na_new+ix,na_new+iz))/dt) *&
((p2(na_new+ix,na_new+iz)-p2prev(na_new+ix,na_new+iz))/dt)

!imagetemp(iz,ix)= &!2.0/(c1(na_new+ix,na_new+iz))*&
 !     p1(na_new+ix,na_new+iz)*p2(na_new+ix,na_new+iz) !gradient result

end forall
image=image+imagetemp

!------------------------------------------------------------------
 ! ! Back propagating boundary values
!------------------------------------------------------------------
if(maxval(p2).gt.10**6) then
write(*,*)'CHECK STABILITY -- VALUES EXPLODING > 10**6'
STOP
end if


!boundary conditions for the tunnel wall and free surface (SH WAVE CASE)
!   u2 --> tauyx
!   w2 --> tauyz
!   p  --> vy
if(tundep.ne.0 .and. wtd.eq.1) then
  w2(1:int(na_new+x_pos-ntungap)+1,int(na_new+ntundep))=0
 w2(1:int(na_new+x_pos-ntungap+1),int(na_new+ntundep+ntunrad))=0
 u2(int(na_new+x_pos-ntungap+1),&
int(na_new+ntundep):int(na_new+ntundep+ntunrad+1))=0
end if
if(free.eq.1 .and. wtd.eq.1) then
w2(1:nx0,na_new+1)=0
end if
!boundary conditions for the free surface (P WAVE CASE)
if(free.eq.1 .and. wtd.eq.2)then
         do ix=1,nx
           p2(ix,na_new+1)=0.0
         enddo
endif

!  applying stored snapshots and boundary value
p2(na_new+1         ,na_new+1:na_new+nz_in)=p2(na_new+1       ,na_new+1:na_new+nz_in)+bordUD(it,:,1)
p2(na_new+nx_in     ,na_new+1:na_new+nz_in)=p2(na_new+nx_in   ,na_new+1:na_new+nz_in)+bordUD(it,:,4)
p2(na_new+2         ,na_new+1:na_new+nz_in)=p2(na_new+2       ,na_new+1:na_new+nz_in)+bordUD(it,:,2)
p2(na_new+nx_in-1   ,na_new+1:na_new+nz_in)=p2(na_new+nx_in-1 ,na_new+1:na_new+nz_in)+bordUD(it,:,3)
if(it.eq.nt-1) then ! snapshot at nt-1
do k=1,nx_in
p2(na_new+k,na_new+1:na_new+nz_in)=p2(na_new+k,na_new+1:na_new+nz_in)+bordUD(1,:,4+k)
end do
end if
p2(na_new+1:na_new+nx_in,      na_new+1)=p2(na_new+1:na_new+nx_in,        na_new+1)+bordRL(it,:,1)
p2(na_new+1:na_new+nx_in,  na_new+nz_in)=p2(na_new+1:na_new+nx_in,    na_new+nz_in)+bordRL(it,:,4)
p2(na_new+1:na_new+nx_in,      na_new+2)=p2(na_new+1:na_new+nx_in,        na_new+2)+bordRL(it,:,2)
p2(na_new+1:na_new+nx_in,na_new+nz_in-1)=p2(na_new+1:na_new+nx_in,  na_new+nz_in-1)+bordRL(it,:,3)
if(it.eq.nt) then ! snapshot at nt
do k=1,nz_in
p2(na_new+1:na_new+nx_in, na_new+k)=p2(na_new+1:na_new+nx_in, na_new+k)+bordRL(1,:,4+k)
end do
end if
!!! adding source
        p2(na_new+x_pos,na_new+z_pos)=p2(na_new+x_pos,na_new+z_pos)-&
        src_sign(it)
p2prev=p2
call move(p2,u2,w2, &
         p2xL, &
         p2xR, &
         p2yL, &
         p2yR, &
         p2xT, &
         p2yT, &
         p2xB, &
         p2yB)
!------------------------------------------------------------------
 ! ! Back propagating data/ data residuals
!------------------------------------------------------------------
if(maxval(p1).gt.10**6) then
write(*,*)'CHECK STABILITY -- VALUES EXPLODING > 10**6'
STOP
end if
!boundary conditions for the tunnel wall and free surface (SH WAVE CASE)
!   u1 --> tauyx
!   w1 --> tauyz
!   p  --> vy

 if(tundep.ne.0 .and. wtd.eq.1) then
   w1(1:int(na_new+x_pos-ntungap)+1,int(na_new+ntundep-1))=0
   w1(1:int(na_new+x_pos-ntungap+1),int(na_new+ntundep+ntunrad+1))=0
   u1(int(na_new+x_pos-ntungap),&
 int(na_new+ntundep-1):int(na_new+ntundep+ntunrad+1))=0
 end if
if(free.eq.1 .and. wtd.eq.1) then
w1(1:nx0,na_new+1)=0
end if
!boundary conditions for the free surface (P WAVE CASE)
if(free.eq.1 .and. wtd.eq.2)then
         do ix=1,nx
           p1(ix,na_new+1)=0.0
         enddo
endif

!!! input data at the receivers

do k=1,ng
p1(na_new+nxg(k,is),na_new+nzg(k,is))= &
    csg(it,k)+p1(na_new+nxg(k,is),na_new+nzg(k,is))
end do

p1prev=p1

!  applying stored snapshots and boundary value
!p1(na_new+1         ,na_new+1:na_new+nz_in)=p1(na_new+1       ,na_new+1:na_new+nz_in)+bordUD(it,:,1)
!p1(na_new+nx_in     ,na_new+1:na_new+nz_in)=p1(na_new+nx_in   ,na_new+1:na_new+nz_in)+bordUD(it,:,4)
!p1(na_new+2         ,na_new+1:na_new+nz_in)=p1(na_new+2       ,na_new+1:na_new+nz_in)+bordUD(it,:,2)
!p1(na_new+nx_in-1   ,na_new+1:na_new+nz_in)=p1(na_new+nx_in-1 ,na_new+1:na_new+nz_in)+bordUD(it,:,3)
if(it.eq.nt-1) then ! snapshot at nt-1
do k=1,nx_in
!p1(na_new+k,na_new+1:na_new+nz_in)=p1(na_new+k,na_new+1:na_new+nz_in)+bordUD(1,:,4+k)
end do
end if
!p1(na_new+1:na_new+nx_in,      na_new+1)=p1(na_new+1:na_new+nx_in,        na_new+1)+bordRL(it,:,1)
!p1(na_new+1:na_new+nx_in,  na_new+nz_in)=p1(na_new+1:na_new+nx_in,    na_new+nz_in)+bordRL(it,:,4)
!p1(na_new+1:na_new+nx_in,      na_new+2)=p1(na_new+1:na_new+nx_in,        na_new+2)+bordRL(it,:,2)
!p1(na_new+1:na_new+nx_in,na_new+nz_in-1)=p1(na_new+1:na_new+nx_in,  na_new+nz_in-1)+bordRL(it,:,3)
if(it.eq.nt) then ! snapshot at nt
do k=1,nz_in
!p1(na_new+1:na_new+nx_in, na_new+k)=p1(na_new+1:na_new+nx_in, na_new+k)+bordRL(1,:,4+k)
end do
end if
call move(p1,u1,w1, &
         p1xL, &
         p1xR, &
         p1yL, &
         p1yR, &
         p1xT, &
         p1yT, &
         p1xB, &
         p1yB)
 
 
!$OMP MASTER
if(it.eq.int(nt/2)) then
call makesu2d(transpose(p1(na_new+1:na_new+nx_in,na_new+1:na_new+nz_in)),&
              trim(out_dir)//'/snap_back',0.0,0.0,dx,dx,nx_in,nz_in,OMP_get_thread_num()+10)
end if
!$OMP END MASTER

end do dt_loop

!  End of Main Loop
!=======================================================
deallocate(c1,c2,cl,den, &
       cl41,cl42,&
       exL,exR, &
       eyT,eyB)

end subroutine back_testing


subroutine move(p2,u1,w1,p2xL,p2xR,p2yL,p2yR,p2xT,p2yT,p2xB,p2yB)
!move/ propagate wavefield in time
implicit none
integer :: ix,iz,ixx,izz
 real,dimension(na_pml,nz0), intent(inout)::p2xL,p2xR, &
       p2yL,p2yR
  real, dimension(nx0,na_pml), intent(inout) :: p2xT,p2yT, &
       p2xB,p2yB
real, dimension(nx0,nz0), intent(inout)::p2,u1,w1


        do iz=nzm10,nzm20
          do ix=nxm10,nxm20
            p2(ix,iz)=p2(ix,iz)+c2(ix,iz)*( &
             c41*( u1(ix  ,iz)-u1(ix-1,iz)+w1(ix,iz  )-w1(ix,iz-1)) &
            +c42*( u1(ix+1,iz)-u1(ix-2,iz)+w1(ix,iz+1)-w1(ix,iz-2)))
          enddo
        enddo


! c p2x

        do iz=na_pml,2,-1
          do ix=nxm10,nxm20
            p2xT(ix,iz)=p2xT(ix,iz)+c2(ix,iz)*( &
            c41*(u1(ix,iz)-u1(ix-1,iz))+c42*(u1(ix+1,iz)-u1(ix-2,iz)))
          enddo
        enddo
        do iz=2,na_pml
          do ix=na_pml,3,-1
            p2xL(ix,iz)=p2xL(ix,iz)*exL(ix,iz)+c2(ix,iz)*( &
            c41*(u1(ix,iz)-u1(ix-1,iz))+c42*(u1(ix+1,iz)-u1(ix-2,iz))) 
          enddo
          ix=2
          p2xL(ix,iz)=p2xL(ix,iz)*exL(ix,iz)+ &
               c2(ix,iz)*(u1(ix,iz)-u1(ix-1,iz))
        enddo
        do ix=nxm20+1,nx-1
          ixx=nx-ix+1
          do iz=2,na_pml
            p2xR(ixx,iz)=p2xR(ixx,iz)*exR(ixx,iz)+c2(ix,iz)*( &
            c41*(u1(ix,iz)-u1(ix-1,iz))+c42*(u1(ix+1,iz)-u1(ix-2,iz)))
          enddo
        enddo
        do iz=nzm20+1,nz-1
           izz=nz-iz+1
	   do ix=nxm10,nxm20
	      p2xB(ix,izz)=p2xB(ix,izz)+c2(ix,iz)*( &
     	   c41*(u1(ix,iz)-u1(ix-1,iz))+c42*(u1(ix+1,iz)-u1(ix-2,iz)))
	   enddo
	enddo
	do ix=na_pml,3,-1
	   do iz=na_pml+1,nz-1
	      p2xL(ix,iz)=p2xL(ix,iz)*exL(ix,iz)+c2(ix,iz)*( &
     	   c41*(u1(ix,iz)-u1(ix-1,iz))+c42*(u1(ix+1,iz)-u1(ix-2,iz)))
	   enddo
	enddo
	ix=2
	do iz=na_pml+1,nz-1
	   p2xL(ix,iz)=p2xL(ix,iz)*exL(ix,iz)+ &
     		c2(ix,iz)*(u1(ix,iz)-u1(ix-1,iz))
	enddo
	do ix=nxm20+1,nx-1
	   ixx=nx-ix+1
	   do iz=na_pml+1,nz-1
	      p2xR(ixx,iz)=p2xR(ixx,iz)*exR(ixx,iz)+c2(ix,iz)*( &
     	   c41*(u1(ix,iz)-u1(ix-1,iz))+c42*(u1(ix+1,iz)-u1(ix-2,iz)))
	   enddo
	enddo
! c
! c p2y
! c
	   do iz=na_pml,3,-1
	   do ix=2,nx-1
	      p2yT(ix,iz)=p2yT(ix,iz)*eyT(ix,iz)+c2(ix,iz)*( &
     	   c41*(w1(ix,iz)-w1(ix,iz-1))+c42*(w1(ix,iz+1)-w1(ix,iz-2)))
	   enddo
	   enddo
	   iz=2
	   do ix=2,nx-1
	   p2yT(ix,iz)=p2yT(ix,iz)*eyT(ix,iz)+ &
     		c2(ix,iz)*(w1(ix,iz)-w1(ix,iz-1))
	   enddo
	do iz=nzm20+1,nz-1
	   izz=nz-iz+1
           do ix=2,nx-1
              p2yB(ix,izz)=p2yB(ix,izz)*eyB(ix,izz)+c2(ix,iz)*( &
          c41*(w1(ix,iz)-w1(ix,iz-1))+c42*(w1(ix,iz+1)-w1(ix,iz-2)))
           enddo
        enddo
	do ix=na_pml,2,-1
	   do iz=nzm10,nzm20
	      p2yL(ix,iz)=p2yL(ix,iz)+c2(ix,iz)*( &
     	   c41*(w1(ix,iz)-w1(ix,iz-1))+c42*(w1(ix,iz+1)-w1(ix,iz-2)))
	   enddo
	enddo
        do ix=nxm20+1,nx-1
	   ixx=nx-ix+1
           do iz=nzm10,nzm20
              p2yR(ixx,iz)=p2yR(ixx,iz)+c2(ix,iz)*( &
          c41*(w1(ix,iz)-w1(ix,iz-1))+c42*(w1(ix,iz+1)-w1(ix,iz-2)))
           enddo
        enddo
! c
	do iz=2,na_pml
	   do ix=2,na_pml
	      p2(ix,iz)=p2xL(ix,iz)+p2yT(ix,iz)
	   enddo
	   do ix=nxm10,nxm20
	      p2(ix,iz)=p2xT(ix,iz)+p2yT(ix,iz)
	   enddo
	   do ix=nxm20+1,nx
	      ixx=nx-ix+1
	      p2(ix,iz)=p2xR(ixx,iz)+p2yT(ix,iz)
	   enddo
	enddo
	do iz=nzm20+1,nz
	   izz=nz-iz+1
           do ix=2,na_pml
              p2(ix,iz)=p2xL(ix,iz)+p2yB(ix,izz)
           enddo
           do ix=nxm10,nxm20
              p2(ix,iz)=p2xB(ix,izz)+p2yB(ix,izz)
           enddo
           do ix=nxm20+1,nx
              ixx=nx-ix+1
              p2(ix,iz)=p2xR(ixx,iz)+p2yB(ix,izz)
           enddo
        enddo   
	do iz=nzm10,nzm20
	   do ix=2,na_pml
	      p2(ix,iz)=p2xL(ix,iz)+p2yL(ix,iz)
	   enddo
	   do ix=nxm20+1,nx
              ixx=nx-ix+1
	      p2(ix,iz)=p2xR(ixx,iz)+p2yR(ixx,iz)
	   enddo
	enddo

     
      do iz=nzm10,nzm20
      do ix=nxm10,nxm20
      u1(ix,iz)=u1(ix,iz)+   &
      cl41(ix,iz)*( p2(ix+1,iz)-p2(ix  ,iz) ) &
      +cl42(ix,iz)*( p2(ix+2,iz)-p2(ix-1,iz) )
      end do
      end do

      do iz=nzm10,nzm20
      do ix=nxm10,nxm20
      w1(ix,iz)=w1(ix,iz)+ &
      cl41(ix,iz)*( p2(ix,iz+1)-p2(ix,iz  ) ) &
     +cl42(ix,iz)*( p2(ix,iz+2)-p2(ix,iz-1) )
      end do
      end do
! ! ! c
! ! c u1
! c
	do iz=na_pml,2,-1
	   do ix=nxm10,nxm20
	      u1(ix,iz)=u1(ix,iz)+ &
      cl41(ix,iz)*( p2(ix+1,iz)-p2(ix  ,iz) ) &
     +cl42(ix,iz)*( p2(ix+2,iz)-p2(ix-1,iz) )  
	   enddo
	enddo
	do iz=nzm20+1,nz-1
	   do ix=nxm10,nxm20
	      u1(ix,iz)=u1(ix,iz)+ &
      cl41(ix,iz)*( p2(ix+1,iz)-p2(ix  ,iz) )&
     +cl42(ix,iz)*( p2(ix+2,iz)-p2(ix-1,iz) )  
	   enddo
	enddo
	do ix=na_pml,2,-1
	   do iz=2,na_pml
	      u1(ix,iz)=u1(ix,iz)*exL(ix,iz)+ &
      cl41(ix,iz)*( p2(ix+1,iz)-p2(ix  ,iz) ) &
     +cl42(ix,iz)*( p2(ix+2,iz)-p2(ix-1,iz) )  
	   enddo
	   do iz=na_pml+1,nz-1
	      u1(ix,iz)=u1(ix,iz)*exL(ix,iz)+ &
      cl41(ix,iz)*( p2(ix+1,iz)-p2(ix  ,iz) ) &
     +cl42(ix,iz)*( p2(ix+2,iz)-p2(ix-1,iz) )  
	   enddo
	enddo
	do ix=nxm20+1,nx-2
	   ixx=nx-ix+1
           do iz=2,na_pml
              u1(ix,iz)=u1(ix,iz)*exR(ixx,iz)+ &
      cl41(ix,iz)*( p2(ix+1,iz)-p2(ix  ,iz) ) &
     +cl42(ix,iz)*( p2(ix+2,iz)-p2(ix-1,iz) )
           enddo
           do iz=na_pml+1,nz-1
              u1(ix,iz)=u1(ix,iz)*exR(ixx,iz)+ &
      cl41(ix,iz)*( p2(ix+1,iz)-p2(ix  ,iz) ) &
     +cl42(ix,iz)*( p2(ix+2,iz)-p2(ix-1,iz) )
           enddo
        enddo
	ix=nx-1
	ixx=nx-ix+1
	do iz=2,na_pml
	   u1(ix,iz)=u1(ix,iz)*exR(ixx,iz)+ &
               cl(ix,iz)*(p2(ix+1,iz)-p2(ix,iz))
	enddo
	do iz=na_pml+1,nz-1
	   u1(ix,iz)=u1(ix,iz)*exR(ixx,iz)+ &
                 cl(ix,iz)*(p2(ix+1,iz)-p2(ix,iz))
	enddo
! c
! c w1
! c
	do iz=na_pml,2,-1
	   do ix=2,nx-1
              w1(ix,iz)=w1(ix,iz)*eyT(ix,iz)+ &
      cl41(ix,iz)*( p2(ix,iz+1)-p2(ix,iz  ) ) &
     +cl42(ix,iz)*( p2(ix,iz+2)-p2(ix,iz-1) )  
	   enddo
	enddo
	do iz=nzm20+1,nz-2
	   izz=nz-iz+1
	   do ix=2,nx-1
	      w1(ix,iz)=w1(ix,iz)*eyB(ix,izz)+ &
      cl41(ix,iz)*( p2(ix,iz+1)-p2(ix,iz  ) ) &
     +cl42(ix,iz)*( p2(ix,iz+2)-p2(ix,iz-1) )
           enddo
        enddo
	iz=nz-1
	izz=nz-iz+1
	do ix=2,nx-1
	   w1(ix,iz)=w1(ix,iz)*eyB(ix,izz)+ &
      cl(ix,iz)*( p2(ix,iz+1)-p2(ix,iz  ) )
	enddo
	do ix=na_pml,2,-1
	   do iz=nzm10,nzm20
	      w1(ix,iz)=w1(ix,iz)+ & 
     cl41(ix,iz)*( p2(ix,iz+1)-p2(ix,iz  ) ) &
     +cl42(ix,iz)*( p2(ix,iz+2)-p2(ix,iz-1) )
           enddo
        enddo
	do ix=nxm20+1,nx
	   do iz=nzm10,nzm20
	      w1(ix,iz)=w1(ix,iz)+ &
     cl41(ix,iz)*( p2(ix,iz+1)-p2(ix,iz  ) ) &
     +cl42(ix,iz)*( p2(ix,iz+2)-p2(ix,iz-1) )
           enddo
	enddo

end subroutine move
!********************************************************************************************

!********************************************************************************************
subroutine do_fwi(f_tol)
! this subroutine does inversion based on global variable inversion_type
! inversion_type	--	'CG' for conjugate gradient
! inversion_type	--	'LBFGS' for LBFGS
! inversion_type	--	'SD' for steepest descent
! in each of the algorithms, the line search is performed by wWolfeLS subroutine
! initial values must be stored in global variable inv_x

! Author: Pawan Bharadwaj
!	  p.b.pisupati@tudelft.nl

! Input:
!	f_tol		-- 	functional tolerence
! 				quits iterations if [(change in functional)/functional].lt.f_tol

! Output:
!	Line search function writes veocity model to disk in each iteration
!	grad_func writes gradient calculated in each iteration

! Global variables used (they must be allocated before calling this subroutine)
!	inv_f,inv_g,inv_p,inv_x,inv_xp


use mpi
use data_global
implicit none
		real, intent(in)		:: f_tol
		real				:: alpha,beta,beta_fr,&
                                                   beta_pr,nmgprevsq,gnorm ! CG parameters

		real, dimension(nx_in*nz_in)	:: y  ! for CG
		real				:: gtp ! for CG
		real,allocatable,dimension(:,:) :: Smat,Ymat ! for lbfgs
		logical				:: go ! when to stop iterations
		real				:: testp ! for lbfgs
		integer				:: M ! storing in lbfgs

!go_fwi
! initializing

itr=0;alpha=0.0;beta=0.0;beta_fr=0.0;beta_pr=0.0;nmgprevsq=0.0;gnorm=0.0;y(:)=0.0;

call grad_func(inv_x,inv_g,inv_f);

if(my_rank.eq.0) then
	write(*,*)''
	write(*,*)'*** inversion flag:',inv_flag
	write(*,*)'*** inversion type:',inversion_type
	write(*,*)'*** functional tolerence:',f_tol
	write(*,*)'frequency band:',iband,'out of',nband_freq,'(',bands_freq(iband,1),'-',bands_freq(iband,2),'Hz)'
	write(*,*)'iteration#, alpha, functional=',itr,alpha,inv_f
endif

call mpi_barrier(mpi_comm_world,ierr)
SELECT CASE (inversion_type)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%
   CASE ('CG') ! conjugate gradient
	
       do itr=1,niteration 
           if(my_rank.eq.0) write(*,*)'===================================================='
           inv_gp=inv_g;inv_fp=inv_f;
           if(itr.gt.1) call test_grad
           gnorm=dot_product(inv_g,inv_g)
           if(itr.eq.1) then
               inv_p(:)=-1.0*inv_g(:);
	       !call update_vel_den(inv_p,alpha) ! first iteration steeest descent 
	       call wWolfeLS(p=inv_p,f0=inv_f,g0=inv_g,x0=inv_x,f_tol=f_tol)
               beta=0;
               if(my_rank.eq.0) write(*,*)'iteration#, alpha, functional=',itr,alpha,inv_f
           else if(itr.gt.1) then
!go_cg     !Calculating beta
           y=inv_g-inv_gp
           nmgprevsq = dot_product(inv_gp,inv_gp)
           beta_pr = dot_product(inv_g,y)/nmgprevsq  ! Polak-Ribiere-Polyak
           beta_fr = dot_product(inv_g,inv_g)/nmgprevsq  ! Fletcher-Reeves  
           if (beta_pr.lt.-1.0*beta_fr) then
               beta = -beta_fr;  
           else if(beta_pr.gt.beta_fr) then
               beta = beta_fr;
           else if(.true.) then 
               beta = beta_pr;
           end if
           if(mod(itr,5).eq.0) beta=0 !CG restart once in every 5 iterations
           inv_p = beta*inv_p - inv_g;
           gtp = dot_product(inv_g,inv_p);
           if(gtp.GE.0) then
	           if(my_rank.eq.0) write(*,*)'not descent direction, quit at iteration',itr;exit
           end if
           inv_xp=inv_x;
	   call wWolfeLS(p=inv_p,f0=inv_f,g0=inv_g,x0=inv_x,f_tol=f_tol)
	   !call update_vel_den(inv_p,alpha)
	   if(((inv_fp-inv_f)/inv_fp).lt.f_tol) then
               if(my_rank.eq.0) write(*,*)'do_fwi: CONVERGED'
           exit 
           endif

	   where(inv_x.gt.vmax) inv_x=vmax
   	   where(inv_x.lt.vmin) inv_x=vmin

           if(my_rank.eq.0) write(*,*)'iteration#, alpha, beta, functional=',itr,alpha,beta,inv_f 
           end if
       end do
!%%%%%%%%%%%%%%%%%%%%%%%%%%%
   CASE ('SD')  !steepest descent
     do itr=1,niteration
           write(*,*)'===================================================='
           inv_gp=inv_g;inv_fp=inv_f;
           inv_p=-1.0*inv_g
           inv_xp=inv_x;
	   call wWolfeLS(p=inv_p,f0=inv_f,g0=inv_g,x0=inv_x,f_tol=f_tol)
           !call update_vel_den(inv_p,alpha)
           if(itr.gt.1) then
           	call test_grad
	        if(((inv_fp-inv_f)/inv_fp).lt.f_tol) then
        		write(*,*)'do_fwi: CONVERGED'
	       		exit	 
	        end if
           endif
           write(*,*)'freq_band#, iteration#, alpha, functional=',iband,itr,alpha,inv_f
      end do 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%
   CASE('LBFGS')
! Simple L-BFGS method 

! input:

!   grad_func - function handle to misfit and gradient calculationsof the form [f,g] = fh(x)
!        where f is the function value, g is the gradient of the same size
!        as the input vector x. 
!   x0 - initial guess
!
!   options.itermax - max iterations [default 10]
!   options.tol     - tolerance on 2-norm of gradient [1e-6]
!   options.M       - history size [5]
!   options.fid     - file id for output [1]
!   options.write   - save iterates to disk [0]

! Author: Pawan Bharadwaj
!	  p.b.pisupati@tudelft.nl

! Date: December, 2012

! You may use this code only under the conditions and terms of the
! license contained in the file LICENSE provided with this source
! code. If you do not agree to these terms you may not use this
! software.


! parameters ---- need to be change

	M         = 5; !(maximum storing of Ymat and Smat to lbfgs method)

allocate(Smat(ninv,M),Ymat(ninv,M))

Smat(:,:) = 0.0;
Ymat(:,:) = 0.0;
go=.true.

! main loop
lbfgs_iterations:do while(go)
    itr=itr+1 
    if(my_rank.eq.0) write(*,*)'===================================================='
    ! compute search direction (p)
    call lbfgs_p(g=-1.0*inv_g,S=Smat,Y=Ymat,p=inv_p)
    ! testing the returned the search direction (reset history if it goes wrong)
    testp = -1.0*dot_product(inv_p,inv_g)/dot_product(inv_g,inv_g);
    
    if(testp.lt.0) then
	if(my_rank.eq.0) write(*,*)'Loss of descent; reset lbfgs history'
        Smat(:,:) = 0.0;
        Ymat(:,:) = 0.0;
	call lbfgs_p(g=-1.0*inv_g,S=Smat,Y=Ymat,p=inv_p)
    endif
    
!go_lbfgs
    inv_gp=inv_g;inv_fp=inv_f;
    inv_xp=inv_x;

    ! linesearch
	
    call wWolfeLS(p=inv_p,f0=inv_f,g0=inv_g,x0=inv_x,f_tol=f_tol,alpha_out=alpha,flag=go)
    
    ! clipping the velocities if blowing

     where(inv_x.lt.vmin) inv_x=vmin
     where(inv_x.gt.vmax) inv_x=vmax

    ! update Smat and Ymat
    Smat=cshift(Smat,shift=1,dim=2)
    Ymat=cshift(Ymat,shift=1,dim=2)
    Smat(:,M)=inv_x-inv_xp
    Ymat(:,M)=inv_g-inv_gp
	      
    write(*,*)'freq_band#, iteration#, alpha, functional=',iband,itr,alpha,inv_f
    
    ! check convergence
    if(itr.gt.niteration) then
	    go=.false.
    endif
enddo lbfgs_iterations
    !call do_lbfgs()
   CASE DEFAULT
      WRITE(*,*)'INVALID INVERSION_TYPE'
      STOP
END SELECT

end subroutine do_fwi

!********************************************************************************************
!********************************************************************************************
subroutine lbfgs_p(g,S,Y,p)
! apply lbfgs inverse Hessian to vector (g to p; gradient to search direction)

! Pawan Bharadwaj
! p.b.pisupati@tudelft.nl

! input:
!    g - vector of length n (gradient)
!    S - history of steps in n x M matrix 
!    Y - history of gradient differences in n x M matrix 

! output
!    p - vector of length n

		integer			:: k, M  !M is storage size for LBFGS
		real,intent(in),dimension(:,:)&
					:: S,Y
		real(dp), intent(in), dimension(:) &
					:: g
		real(dp),intent(out)	:: p(:)
		real			:: alpha(size(S,2)),rho(size(S,2)),q(size(g))
		real			:: a,beta,temp

if(size(S,1).ne.size(Y,1) .or. size(S,1) .ne. size(g)) then
	write(*,*)'lbfgs_p: size of S and Y (or) g not the same'; stop
endif
if(size(g).ne.size(p)) then
	write(*,*)'lbfgs_p: size of g and p should be same'; stop
endif
if(size(S,2).ne.size(Y,2)) then
	write(*,*)'lbfgs_p: Y and S should have same length'; stop
endif
!go_lbfgs_p
	p(:)=0.0;
	M = size(S,2);
	alpha(:)=0.0; rho(:)=0.0;
	if(maxval(S)*maxval(Y).ne.0) then
		do k = 1,M
			temp=dot_product(Y(:,k),S(:,k));
			if(temp.ne.0) rho(k) = 1.0/temp
		enddo
	endif
	q = g;
	! first recursion
	if(maxval(S)*maxval(Y).ne.0) then
		do k = M,1,-1
			alpha(k) = rho(k)*dot_product(S(:,k),q);
			q        = q - alpha(k)*Y(:,k);
		enddo
	endif
	! apply `initial' Hessian
	if(M.gt.0 .and. maxval(S)*maxval(Y).ne.0) then
		temp=dot_product(Y(:,M),Y(:,M));
		if(temp.eq.0) then 
			write(*,*)'lbfgs_p: division by zero'	
			stop
		endif
		a = (dot_product(Y(:,M),S(:,M))/temp);
	else
		a = 1.0/sum(abs(g));
	endif
	p = a*q;
	! second recursion
	if(maxval(S)*maxval(Y).ne.0) then
		do k = 1, M
			beta = rho(k)*(dot_product(Y(:,k),p));
			p    = p + (alpha(k) - beta)*S(:,k);
		enddo
	endif
write(*,*)'maxval(p)',maxval(p)

end subroutine lbfgs_p

!********************************************************************************************
!********************************************************************************************

subroutine grad_func(x,gradient,functional)
! this subroutine returns gradient and functional value at x.
! which gradient of which functional ? determined by a global variable called -- inv_flag

! Author : Pawan Bharadwaj
!	   p.b.pisupati@tudelft.nl

! Inputs:
!		x	-->	inversion variable

! Outputs:
!		f	--> 	functional at x (optional)
!		g	-->	gradient at x (optional)

real(dp), intent(in), optional		:: x(ninv) ! inversion variable
real(dp), intent(out), optional		:: gradient(ninv) 
real(dp), intent(out), optional		:: functional

real 					:: rtm(nz_in,nx_in),f,ctemp1(nz_in*nx_in),ctemp2(nz_in*nx_in)&
						,g2(nz_in*nx_in),f2
character(len=200)			:: junk
real(dp),allocatable,dimension(:)	:: functionals
real, dimension(:,:,:),allocatable	:: csgs_obs_AllFreq

if(present(gradient)) gradient(:)=0.0
if(present(functional)) functional=0.0
csgs_cal(:,:,:)=0.0

select case(abs(inv_flag))

!%%%%%%%%%%%%%%%%%%%%%%%%%%%
! REVERSE TIME MIGRATION
! in this case gradient is the RTM image formed by 
! which is given by sum(all shots of all stages) 
! x and f makes no sense
case(1)
	allocate(illu(nz_in,nx_in,nstage),illu_rec(nz_in,nx_in),illu_src(nz_in,nx_in))
	do istg=1,nstage/nump
		istage=1+my_rank+(istg-1)*ntsage/nump
		write(*,*)'stage#	',istage
		call do_rtm(rtm,0)
		gradient(1:nx_in*nz_in)=gradient(1:nx_in*nz_in)+pack(rtm,.true.)
	enddo
!go_grad_func
!%%%%%%%%%%%%%%%%%%%%%%%%%%%
! CONVENTIONAL FWI
! gradient and functional of only iband
! functional :  least square data residual
! gradient   :  sum over all individual source gradients of all stages 
case(2)
        c=sngl(x)
	!call put_tunnel(model_out=c,model_in=sngl(x),value=0.0)

	if(present(gradient)) then
		allocate(illu(nz_in,nx_in,nstage),illu_rec(nz_in,nx_in),illu_src(nz_in,nx_in))
		do istg=1,nstage/nump
			istage=1+my_rank+(istg-1)*ntsage/nump
			write(*,*)'stage#	',istage
			call do_rtm(rtm,1)
			gradient=pack(rtm,.true.)
		enddo
		call mpi_add_1ddp(gradient)
		call mpi_add_3d(csgs_cal)
		call filename(junk,trim(out_dir)//'/grad_band#',iband,'_itr#',itr)
		call makesu2d(reshape(sngl(gradient),(/nz_in,nx_in/)),junk,0.0,0.0,dx,dx,nx_in,nz_in,10)

	call writebin3d(trim(out_dir)//'/csgs_test1.bin',csgs_cal,nt,ng,ns) !saving recorded data
	endif

	if(present(functional)) then
		if(present(gradient).eqv..false.) then
			call do_fw_modelling(0)
		endif
		do istg=1,nstage
			functional=functional+distance3d(csgs_cal(:,:,ns/nstage*(istg-1)+1:ns/nstage*istg),&
        	            csgs_obs(:,:,ns/nstage*(istg-1)+1:ns/nstage*istg));
		enddo
	call writebin3d(trim(out_dir)//'/csgs_test2.bin',csgs_cal,nt,ng,ns) !saving recorded data
	endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%
! gradient and functional of only istage
! functional :  returns least square functional of a particular istage
! gradient   :  returns gradient of a particular stage (istage)
case(4)
	
	call put_tunnel(model_out=c,model_in=sngl(x),value=0.0)
	if(present(gradient)) then
		allocate(illu(nz_in,nx_in,nstage),illu_rec(nz_in,nx_in),illu_src(nz_in,nx_in))
		write(*,*)'stage#	',istage
		call do_rtm(rtm,1)
       		gradient=gradient+pack(rtm,.true.)
	endif
	if(present(functional)) then
		if(present(gradient).eqv..false.) then
			call do_fw_modelling(0)
		endif
		f= distance3d(csgs_cal(:,:,ns/nstage*(istage-1)+1:ns/nstage*istage),&
                    csgs_obs(:,:,ns/nstage*(istage-1)+1:ns/nstage*istage));
		functional=f
	endif

if(my_rank.eq.0) then
	call filename(junk,trim(out_dir)//'/csgs_cal_band#',iband,'_itr#',itr)
	call writebin3d(junk,csgs_cal,nt,ng,ns) !saving recorded data
endif
call mpi_add_3d(csgs_cal)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%
!go_grad_func
! Wim's idea of multi dimensional problem
! 
! J(stage) = Jls(stage) + 1/2(c(stage)-c0)^2
! functional = sum(J(stage))
case(6)
	c=sngl(x(nx_in*nz_in+1:nx_in*nz_in*(nstage+1)))
	if(present(gradient)) then
	allocate(illu(nz_in,nx_in,nstage),illu_rec(nz_in,nx_in),illu_src(nz_in,nx_in))
		do istg=1,nstage/nump
			istage=1+my_rank+(istg-1)*ntsage/nump
			ctemp1=x(1:nz_in*nx_in)
			ctemp2=x(nx_in*nz_in*(istage)+1:nx_in*nz_in*(istage+1))
			write(*,*)'stage#	',istage
			call do_rtm(rtm,1)
			gradient(nx_in*nz_in*istage+1:nx_in*nz_in*(istage+1))=pack(rtm,.true.)+&
                                        eps*(ctemp2-ctemp1)
		end do
		call mpi_add_1ddp(gradient)
		call mpi_add_3d(csgs_cal)
		if(my_rank.eq.0) then
			do istg=1,nstage
				call filename(junk,trim(out_dir)//'/grad_band#',iband,'_stage#',istg,'_itr#',itr)
				call makesu2d(sngl(gradient&
				    (nx_in*nz_in*istg+1:nx_in*nz_in*(istg+1))),junk,0.0,0.0,dx,dx,nx_in,nz_in,10)
			enddo
		endif

		do istg=1,nstage
			ctemp1=x(1:nz_in*nx_in)
			ctemp2=x(nx_in*nz_in*(istg)+1:nx_in*nz_in*(istg+1))
                	gradient(1:nx_in*nz_in)=gradient(1:nx_in*nz_in)+eps*(ctemp1-ctemp2)
		enddo
		if(my_rank.eq.0) then
			call filename(junk,trim(out_dir)//'/grad_bkmodel_itr#',itr,'')
			call makesu2d(sngl(gradient),junk,0.0,0.0,dx,dx,nx_in,nz_in,10)
		endif

	endif

	if(present(functional)) then
		if(present(gradient).eqv..false.) then
			call do_fw_modelling(0)
		endif
		do istg=1,nstage
			ctemp1=x(1:nz_in*nx_in)
			ctemp2=x(nx_in*nz_in*(istg)+1:nx_in*nz_in*(istg+1))
			f=f+1.0/2.0*(distance3d(csgs_cal(:,:,ns/nstage*(istg-1)+1:ns/nstage*istg),&
                            csgs_obs(:,:,ns/nstage*(istg-1)+1:ns/nstage*istg)) &
                            +eps*dot_product((ctemp2-ctemp1),(ctemp2-ctemp1)))
		enddo
		functional=f;
	endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Multi dimensional problem - \sum{stage} [ LS functional of stage ]  + &
!		1/2 * \sum{stage1} \sum{stage2} ||mstage1-mstage2||2  

case(7)
	c=sngl(x)
	if(present(gradient)) then
	allocate(illu(nz_in,nx_in,nstage),illu_rec(nz_in,nx_in),illu_src(nz_in,nx_in))
		do istg=1,nstage/nump
			istage=1+my_rank+(istg-1)*ntsage/nump
			g2(:)=0.0;
 			do k=1,nstage
				ctemp1=x(nx_in*nz_in*(istage-1)+1:nx_in*nz_in*(istage))
				ctemp2=x(nx_in*nz_in*(k-1)+1:nx_in*nz_in*(k))
				call fg_model(m0=ctemp2,mi=ctemp1,g=g2,flag=1)
				
			enddo
			write(*,*)'stage#	',istage
			call do_rtm(rtm,1)
			gradient(nx_in*nz_in*(istage-1)+1:nx_in*nz_in*(istage))=pack(rtm,.true.)+&
                                        eps*g2
		end do
		call mpi_add_1ddp(gradient)
		call mpi_add_3d(csgs_cal)
		if(my_rank.eq.0) then
			do istg=1,nstage
				call filename(junk,trim(out_dir)//'/grad_band#',iband,'_stage#',istg,'_itr#',itr)
				call makesu2d(sngl(gradient&
				    (nx_in*nz_in*(istg-1)+1:nx_in*nz_in*(istg))),junk,0.0,0.0,dx,dx,nx_in,nz_in,10)
			enddo

		endif

	endif

	if(present(functional)) then
		if(present(gradient).eqv..false.) then
			call do_fw_modelling(0)
		endif
		do istg=1,nstage
			f2=0.0;
			do k=1,nstage
				ctemp1=x(nx_in*nz_in*(istg-1)+1:nx_in*nz_in*(istg))
				ctemp2=x(nx_in*nz_in*(k-1)+1:nx_in*nz_in*(k))
				call fg_model(m0=ctemp2,mi=ctemp1,f=f2,flag=1)
			enddo
			f=f+1.0/2.0*(distance3d(csgs_cal(:,:,ns/nstage*(istg-1)+1:ns/nstage*istg),&
                            csgs_obs(:,:,ns/nstage*(istg-1)+1:ns/nstage*istg)))
		enddo
		functional=f+eps*f2;
	endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%
case default
	write(*,*)'GRAD_FUNC: INVALID INVERSION_FLAG'
	stop
end select
if(present(gradient)) then
	deallocate(illu,illu_rec,illu_src)
end if

end subroutine grad_func
!********************************************************************************************

!********************************************************************************************
subroutine fg_model(m0,mi,g,f,flag)
! this subroutine can calculate functional and gradient of different 
! cost functions which measure similarity between two models
! 
! input:
! 		flag 		: 1 - ||m0-mi||2
!				: 2 - need to be added
!		m0		: reference model
!		mi		: variable model
! m0 and mi should have same dimensions
! output (g and f are modified according to )
! 		g		: g + gradient of f wrt mi (optional)
!		f		: f + functional value (optional)

! if m0 or mi have zeros then gradients of those points are not calculated
! and also not included while calculating functional

implicit none

	real, intent(in)		:: m0(:),mi(:)
	real, intent(inout), optional	:: f,g(size(m0))
	integer, intent(in)		:: flag
	integer				:: n,i

n=size(m0)
if(n.ne.size(mi)) then
	write(*,*)'fg_model: reference and variable model should have same dimension'
	stop	
endif

if(present(f)) then
	do i=1,n
		if(mi(i).ne.0.0 .and. m0(i).ne.0.0) then
			f=f+1.0/2.0*(mi(i)-m0(i))**2
		endif
	enddo
endif

if(present(g)) then
	do i=1,n
		if(mi(i).ne.0.0 .and. m0(i).ne.0.0) then
			g(i)=mi(i)-m0(i)
		endif
	enddo
endif


end subroutine fg_model
!********************************************************************************************

!********************************************************************************************
subroutine test_grad
implicit none
integer ix,iz
real test(ninv),test2 
!gradient test, Wim's idea .. testing if (Fn+1-Fn)/2(gn+1+gn)(mn+1-mn) is close to 1
do ix=1,nx_in*nz_in
      test(ix)=(inv_g(ix)+inv_gp(ix))/2.0*(inv_x(ix)-inv_xp(ix))
end do
!go_test_grad
test2=sum(test)/(inv_f-inv_fp)
!call makesu2d(unpack(test,TRUEMAT,ZEROMAT),'./output/test_grad',0.0,0.0,dx,dx,nx_in,nz_in,10)
write(*,*)'Testing calculated gradient (value should be close to 1)', test2
end subroutine test_grad
!********************************************************************************************
!********************************************************************************************
subroutine do_rtm(rtm, method)
! this subroutine calculates the gradient for all shots in a particular stage
! method	--0 back propagates obsserved data (csgs_obs) -- for rtm
! method	--1 back propagates csg_obs-csg_calculated 
!		    and gradients for different source positions 
!		    simply added (for conventional LS FWI)
!                   also gradients will be multiplied with inverse of approximate Hessian (illu)
! method	--2 back propagates csg_obs-csg_calculated 
!		    and gradients for different source positions 
!		    added after weighting with respective LS errors (for inversion_flag 3)
!                   also gradients will be multiplied with inverse of approximate Hessian (illu)

use data_global
use omp_lib
implicit none
integer				:: is,i,j,it,ix,iz,ig
integer , intent(in)		:: method
real,intent(out)		:: rtm(nz_in,nx_in)
real				:: csg_impulse(nt,ng),src_temp(nt)
real				:: bord_UD(nt,nz_in,4+nx_in),bord_RL(nt,nx_in,4+nz_in)
real, dimension(:,:,:,:), allocatable &
				:: snaps
real,dimension(nz_in,nx_in)	:: rtm_temp,illu_src_temp
real				:: residual(nt,ng) 
character(len=200)		:: junk
illu(:,:,istage)=0
illu_src(:,:)=0
illu_rec(:,:)=0
do it=1,nt
	csg_impulse(nt-it+1,:)=src_sign(it)        
enddo

rtm(:,:)=0
!applying bounds to the velocity values
where(c.gt.vmax) c=vmax
where(c.lt.vmin.and.c.ne.0) c=vmin
allocate(snaps(nt,nx_in,nz_in,int(ns/nstage)))
if(illu_flag.eq.1) then
	write(*,*)'calculating receiver illumination ...' 
	!$OMP PARALLEL NUM_THREADS(num_threads) DEFAULT(SHARED) &
	!$OMP PRIVATE(it,ix,iz,src_temp)
	!$OMP DO
	do ig=1,ng
	    if(method.gt.0) then
	      call do_bw_filter1d(src_temp,src_sign,bands_freq(iband,1),bands_freq(iband,2),dt)
	      call forw(vel_flag=0,is=ns/nstage*(istage-1)+1,illu=illu_rec,ig=ig,source_wav=src_temp)
	    else if(method.eq.0) then
	      call forw(vel_flag=-1,is=ns/nstage*(istage-1)+1,illu=illu_rec,ig=ig,source_wav=src_sign)
	   endif
	enddo	
	!$OMP END DO 
	!$OMP END PARALLEL
endif
write(*,*)'calculating gradients for different sources ...'
csgs_cal(:,:,:)=0.0;
!$OMP PARALLEL NUM_THREADS(num_threads) DEFAULT(SHARED) &
!$OMP PRIVATE(it,junk,rtm_temp,ix,iz,illu_src_temp,src_temp)
!$OMP DO
do is=ns/nstage*(istage-1)+1,ns/nstage*istage
    illu_src_temp(:,:)=0
    rtm_temp(:,:)=0
    !write(*,*)'source position(x,z)=',xs(is),'m,',zs(is),'m'
    if(method.gt.0) then
        if(iband.ne.0) then
        call do_bw_filter1d(src_temp,src_sign,bands_freq(iband,1),bands_freq(iband,2),dt)
        endif
        call forw(vel_flag=0,csg=csgs_cal(:,:,is),illu=illu_src_temp,&
        is=is,source_wav=src_temp,snaps=snaps(:,:,:,is-ns/nstage*(istage-1)))
        call back(vel_flag=0,grad=rtm_temp,snaps=snaps(:,:,:,is-ns/nstage*(istage-1))&
                    ,is1=is,csg1=csgs_cal(:,:,is)-csgs_obs(:,:,is))
    else if(method.eq.0) then
        call forw(vel_flag=-1,csg=csgs_cal(:,:,is),illu=illu_src_temp,is=is,source_wav=src_sign,&
                         snaps=snaps(:,:,:,is-ns/nstage*(istage-1)))
        call back(vel_flag=-1,grad=rtm_temp,snaps=snaps(:,:,:,is-ns/nstage*(istage-1)),&
                       is1=is,csg1=csgs_obs(:,:,is))
    end if
    
! illumination compensation
! call illu_muting(rtm_temp,illu_src_temp*illu_rec,0.01);
    do ix=1,nx_in
    do iz=1,nz_in
        if(illu_src_temp(iz,ix).gt.0) rtm_temp(iz,ix)=rtm_temp(iz,ix)/illu_src_temp(iz,ix)
        if(illu_rec(iz,ix).gt.0) rtm_temp(iz,ix)=rtm_temp(iz,ix)/illu_rec(iz,ix)
    enddo
    enddo
    call mute_grad(rtm_temp,1,is)
 !   call filename(junk,trim(out_dir)//'/rtm_shot#',is,'')
 !   call makesu2d(rtm_temp,junk,0.0,0.0,dx,dx,nx_in,nz_in,OMP_get_thread_num()+10)
    if(method.le.1) then
        rtm=rtm+rtm_temp
    else if(method.eq.2) then !gradient of method 3 functional
        rtm=rtm+rtm_temp/(distance2d(csgs_cal(:,:,is),csgs_obs(:,:,is)))
    endif
    illu_src=illu_src+illu_src_temp
end do
!go_rtm
!$OMP END DO 
!$OMP END PARALLEL
deallocate(snaps)
illu(:,:,istage)=illu_src*illu_rec
!call makesu2d(illu(:,:,istage),trim(out_dir)//'/illu',0.0,0.0,dx,dx,nx_in,nz_in,10)
!call makesu2d(illu_src,trim(out_dir)//'/illu_src',0.0,0.0,dx,dx,nx_in,nz_in,10)
!call makesu2d(illu_rec,trim(out_dir)//'/illu_rec',0.0,0.0,dx,dx,nx_in,nz_in,10)
!call filename(junk,trim(out_dir)//'/csgs_cal_band#',iband,'_itr#',itr)
!call writebin3d(junk,csgs_cal,nt,ng,ns) !saving recorded data
call filename(junk,trim(out_dir)//'/csgs_cal_band#',iband,'_itr#',itr)
!call writebin3d(junk,csgs_cal,nt,ng,ns) !saving observed data for this band
call filename(junk,trim(out_dir)//'/residual_band#',iband,'_itr#',itr)
!call writebin3d(junk,csgs_cal-csgs_obs,nt,ng,ns) !saving residual
call filename(junk,trim(out_dir)//'/grad_allshot_stage#',istage,'_itr#',itr)
!call makesu2d(rtm,junk,0.0,0.0,dx,dx,nx_in,nz_in,10)
if(method.eq.0) then
write(*,*)'===================================================='
end if
end subroutine do_rtm
!********************************************************************************************

!********************************************************************************************
subroutine do_rtm2(rtm)
use data_global
use omp_lib
implicit none
integer :: is,i,j,it,ix,iz,ig
real,intent(out) :: rtm(nz_in,nx_in)
real :: bord_UD(nt,nz_in,4+nx_in,ns),bord_RL(nt,nx_in,4+nz_in,ns)
real,dimension(nz_in,nx_in) :: rtm_temp
! illu_src	::	source illumination
! illu_rec	::	receiver illumination
character(len=200) :: junk
illu_src(:,:)=0
rtm(:,:)=0
!$OMP PARALLEL NUM_THREADS(num_threads) DEFAULT(SHARED) &
!$OMP PRIVATE(it,junk,rtm_temp,ix,iz)
!$OMP DO
!go_rtm2
do is=1,ns
        rtm_temp(:,:)=0
        write(*,*)'source position(x,z)=',xs(is),'m,',zs(is),'m'
        call forw(vel_flag=-1,is=is,illu=illu_src,&
             source_wav=src_sign,bord_UD=bord_UD(:,:,:,is),bord_RL=bord_RL(:,:,:,is))
        call back_testing(vel_flag=-1,image=rtm_temp,&
             bordUD=bord_UD(:,:,:,is),bordRL=bord_RL(:,:,:,is),is=is,csg=csgs_obs(:,:,is))
        call filename(junk,trim(out_dir)//'/rtm_shot#',is,'')
        call makesu2d(rtm_temp,junk,0.0,0.0,dx,dx,nx_in,nz_in,OMP_get_thread_num()+10)
        !call writebin2d(junk,rtm_temp,nz_in,nx_in,OMP_get_thread_num()+10)
        rtm=rtm+rtm_temp
end do
!$OMP END DO 
!$OMP END PARALLEL

!$OMP PARALLEL NUM_THREADS(num_threads) DEFAULT(SHARED) &
!$OMP PRIVATE(it,ix,iz)
!$OMP DO
do ig=1,ng
       call forw(vel_flag=-1,is=1,illu=illu_rec,ig=ig,source_wav=src_sign)
enddo
!$OMP END DO 
!$OMP END PARALLEL
rtm=1.0/(illu_src+illu_rec)*rtm !illumination compensation
call filename(junk,trim(out_dir)//'/grad_allshot_itr#',itr,'')
call makesu2d(illu_src,trim(out_dir)//'/illu_src',0.0,0.0,dx,dx,nx_in,nz_in,10)
call makesu2d(illu_rec,trim(out_dir)//'/illu_rec',0.0,0.0,dx,dx,nx_in,nz_in,10)
call makesu2d(rtm,junk,0.0,0.0,dx,dx,nx_in,nz_in,10)
end subroutine do_rtm2
!********************************************************************************************

!****************************************************************************************
subroutine update_vel_den(pin,alpha)
! This subroutine updates inv_x which is global variable for inversion
! if p is present the step_length is calculated by parabolic approximation

           implicit none
           real               :: kappa 
           real(dp), intent(in) &
                              :: pin(ninv)
           real,intent(out)   :: alpha
           real               :: csgs_temp(nt,ng,ns)
           integer            :: is,istg,nxnz
           real               :: mu1,res1,mu2,res2
           character(len=200) :: junk
           optional           :: pin, alpha
!pin -- in which direction  should i move
!calculation of step length using  
!Pica, A., J. P. Diet, and A. , 1990, Nonlinear inversion of seismic re-
!flection data in a laterally invariant medium: Geophysics, 55, 284–292
if(present(pin)) then
select case(abs(inv_flag))
!########################
   case default   
   kappa=(vmax)/(maxval(abs(pin))*10)
   write(*,*)'calculating step length..., kappa=',kappa
   !go_alpha
   call put_tunnel(model_out=c,model_in=sngl(inv_x+kappa*pin),value=tunvel)
   !$OMP PARALLEL NUM_THREADS(num_threads) DEFAULT(SHARED) &
   !$OMP PRIVATE(junk)
   !$OMP DO
   do is=1,ns
      call forw(vel_flag=0,csg=csgs_temp(:,:,is),is=is,source_wav=src_sign)
      if(iband.ne.0) then
           call do_bw_filter(csgs_temp(:,:,is),csgs_temp(:,:,is),bands_freq(iband,1),bands_freq(iband,2),dt)
      endif
   end do
   !$OMP END DO 
   !$OMP END PARALLEL
   res1=dot_product(pack(csgs_temp-csgs_cal,.true.),pack(csgs_temp-csgs_cal,.true.))
   res1=res1/kappa/kappa
   mu1=dot_product(pack(csgs_cal-csgs_obs,.true.),pack(csgs_cal-csgs_temp,.true.)/kappa)
   mu2=0;res2=0;
   ! mu=dot_product(pin,pin) !testing
   alpha=(mu1+eps*mu2)/(res1+eps*res2)
   call put_tunnel(model_out=c,model_in=sngl(inv_x),value=tunvel)
   inv_x=inv_x+alpha*pin
   write(*,*)'mu, res',mu1,res1
   call filename(junk,trim(out_dir)//'/vel_band#',iband,'_itr#',itr)
   call makesu2d(reshape(sngl(inv_x),(/nz_in,nx_in/))&
                  ,junk,0.0,0.0,dx,dx,nx_in,nz_in,10)
!########################
   case(4)   
   kappa=(vmax)/(maxval(abs(pin))*100)
   write(*,*)'calculating step length..., kappa=',kappa
   !go_alpha
   c(1+(istage-1)*nx_in*nz_in:nx_in*nz_in*istage)=inv_x+kappa*pin
   csgs_temp(:,:,:)=0.0;
   !$OMP PARALLEL NUM_THREADS(num_threads) DEFAULT(SHARED) &
   !$OMP PRIVATE(junk)
   !$OMP DO
   do is=ns/nstage*(istage-1)+1,ns/nstage*istage
      call forw(vel_flag=0,csg=csgs_temp(:,:,is),is=is,source_wav=src_sign)
      if(iband.ne.0) then
           call do_bw_filter(csgs_temp(:,:,is),csgs_temp(:,:,is),bands_freq(iband,1),bands_freq(iband,2),dt)
      endif
   end do
   !$OMP END DO 
   !$OMP END PARALLEL
   call mpi_add_3d(csgs_temp)
   res1=dot_product(pack(csgs_temp(:,:,ns/nstage*(istage-1)+1:ns/nstage*istage)-&
   csgs_cal(:,:,ns/nstage*(istage-1)+1:ns/nstage*istage),.true.),&
   pack(csgs_temp(:,:,ns/nstage*(istage-1)+1:ns/nstage*istage)-&
   csgs_cal(:,:,ns/nstage*(istage-1)+1:ns/nstage*istage),.true.))
   res1=res1/kappa/kappa
   mu1=dot_product(pack(csgs_cal(:,:,ns/nstage*(istage-1)+1:ns/nstage*istage)-&
   csgs_obs(:,:,ns/nstage*(istage-1)+1:ns/nstage*istage),.true.),&
   pack(csgs_cal(:,:,ns/nstage*(istage-1)+1:ns/nstage*istage)-&
   csgs_temp(:,:,ns/nstage*(istage-1)+1:ns/nstage*istage),.true.)/kappa)

   mu2=dot_product(inv_x-bkmodel,pin);
   res2=dot_product(pin,pin)
   ! mu=dot_product(pin,pin) !testing
   alpha=(mu1+eps*mu2)/(res1+eps*res2)
   c(1+(istage-1)*nx_in*nz_in:nx_in*nz_in*istage)=inv_x
   inv_x=inv_x+alpha*pin
   write(*,*)'mu, res, alpha',mu1,res1,alpha
   !call filename(junk,trim(out_dir)//'/vel_band#',iband,'_itr#',itr)
   !call makesu2d(reshape(sngl(inv_x),(/nz_in,nx_in/))&
   !               ,junk,0.0,0.0,dx,dx,nx_in,nz_in,10)
!########################

   case(6)
   nxnz=nx_in*nz_in
   kappa=(vmax)/(maxval(abs(pin(nxnz+1:nxnz*nstage)))*10)
   write(*,*)'calculating step length..., kappa=',kappa
   !go_alpha
   c=sngl(inv_x(nxnz+1:nxnz*nstage)+kappa*pin(nxnz+1:nxnz*nstage))
   !$OMP PARALLEL NUM_THREADS(num_threads) DEFAULT(SHARED) &
   !$OMP PRIVATE(junk)
   !$OMP DO
   do is=1,ns
      call forw(vel_flag=0,csg=csgs_temp(:,:,is),is=is,source_wav=src_sign)
      if(iband.ne.0) then
           call do_bw_filter(csgs_temp(:,:,is),csgs_temp(:,:,is),bands_freq(iband,1),bands_freq(iband,2),dt)
      endif
   end do
   !$OMP END DO 
   !$OMP END PARALLEL
   res1=dot_product(pack(csgs_temp-csgs_cal,.true.),pack(csgs_temp-csgs_cal,.true.))
   res1=res1/kappa/kappa
   mu1=dot_product(pack(csgs_cal-csgs_obs,.true.),pack(csgs_cal-csgs_temp,.true.)/kappa)
   mu2=0;res2=0;
   do istg=1,nstage
        mu2=mu2+dot_product((inv_x(nxnz*istg+1:nxnz*(istg+1))-inv_x(1:nxnz)),&
                    (pin(nxnz*istg+1:nxnz*(istg+1))-pin(1:nxnz)))
        res2=res2+dot_product((pin(nxnz*istg+1:nxnz*(istg+1))-pin(1:nxnz)),&
                    (pin(nxnz*istg+1:nxnz*(istg+1))-pin(1:nxnz)))
   enddo
   alpha=(mu1-eps*mu2)/(res1+eps*res2)
   inv_x=inv_x+alpha*pin
   write(*,*)'mu1, res1, mu2, res2',mu1,res1,mu2,res2
   do istg=0,nstage
        call filename(junk,trim(out_dir)//'/vel_stage#',istg,'_itr#',itr)
        call makesu2d(reshape(sngl(inv_x(nxnz*(istg)+1:nxnz*(istg+1))),(/nz_in,nx_in/))&
                        ,junk,0.0,0.0,dx,dx,nx_in,nz_in,10)
        call filename(junk,trim(out_dir)//'/grad_stage#',istg,'_itr#',itr)
        call makesu2d(reshape(sngl(inv_g(nxnz*(istg)+1:nxnz*(istg+1))),(/nz_in,nx_in/))&
                        ,junk,0.0,0.0,dx,dx,nx_in,nz_in,10)
   enddo
end select
endif
end subroutine update_vel_den
!********************************************************************************************
!********************************************************************************************
subroutine mute_grad(g,mode,is)
! this subroutine preconditions the calculated gradient 
! input  -- g -- gradient  (nx_in,nz_in)
! mode 	 	-- 1 mute around sources
!		-- 2 mute around receivers 
!		-- 3 both around sources and receivers ..
!		-- 4 ??
! is		-- around which source and receivers of which source.. 
real, intent(inout)		:: g(nz_in,nx_in);
real				:: grad(nz_in,nx_in)
integer, intent(in)		:: is,mode
integer				:: i,ix,iz
! radius of the circle to be muted .. not in meter.. but in #grid point
real				:: r
r=5
do ix=1,nx_in
do iz=1,nz_in
if((((ix-nxs(is))**2)+(iz-nzs(is))**2).lt.r**2) then
g(iz,ix)=0;
endif
!do ig=1,ng
!if((((ix-nxg(ig,is))**2)+(iz-nzg(ig,is))**2).lt.r**2) then
!g(iz,ix)=0
!endif
!enddo
enddo
enddo
end subroutine mute_grad
!********************************************************************************************
!********************************************************************************************
subroutine put_tunnel(model_out,model_in,value)
! this subroutine inserts a tunnel in model for a given stage..
! model is a vector of nx_in*nz_in*nstage
! value -- value to be put in the tunnel
         real,intent(inout), dimension(nz_in*nx_in*nstage)&
                                          :: model_out
         real, intent(in), dimension(nz_in*nx_in) &
                                          :: model_in
         real,intent(in)                  :: value
         real, dimension(nz_in,nx_in)     :: temp
         integer                          :: ix,iz,istg
do istg=1,nstage
temp=reshape(model_in,(/nz_in,nx_in/))
    if(tundep.ne.0) then
        ntundep=int(tundep/dx)
        ntunrad=int(tunrad/dx)
        ntungap=int(tungap/dx)
        do ix=1,nxs(1+ns/nstage*(istg-1))-ntungap
                do iz=ntundep,ntundep+ntunrad
                        temp(iz,ix)=value
                end do
        end do
    endif
model_out(nx_in*nz_in*(istg-1)+1:nx_in*nz_in*(istg))=pack(temp,.true.)
enddo
end subroutine
!********************************************************************************************
!********************************************************************************************

subroutine wWolfeLS(x0,f0,g0,p,f_tol,alpha_out,flag)
! Simple Wolfe linesearch, adapted from
! (http://cs.nyu.edu/overton/mstheses/skajaa/msthesis.pdf, algorihtm 3).

! Authr : Pawan Bharadwaj
!	  p.b.pisupati@tudelft.nl

! Function used for calculation of f and g -- > grad_func

! Inputs:
!		x0 --> initial value of variable
!		f0 --> initial functional
!		g0 --> gradient at x0
!		p  --> which direction to move ?
!		f_tol --> functional tolerence when quit line search.

! Outputs: (as global variables)
!		inv_f  --> functional value after line search is complete
!		inv_g  --> gradient at inv_x
!		inv_x  --> x0+alpha*p
!		alpha  --> step length optional
! 		flag   --> flag =.true. if line search is success
!			        =.false. if line search fails 
!		
use data_global
use mpi
implicit none
	integer			::	lsiter,istg ! number of functional evaluations
	!real			::	c1,c2
	logical, intent(inout), optional&
				::      flag
	logical			::	go
	real(dp)		::	alpha,mu
	real			::	f_tol
	real(dp), intent(in)	::	f0,g0(ninv),p(ninv),x0(ninv)
	real(dp)		::	ft,gt(ninv),fp
	character(len=200)	::	junk
	real, intent(out),optional&
				::	alpha_out
	
! go_line
go=.true.
lsiter = 0;
mu = 0;
alpha=vmin/10/maxval(p); ! initial guess
ft=f0;fp=f0;

do while(go)

    if(lsiter.lt.40) then
        call grad_func(x=x0+alpha*p,functional=ft)
        lsiter = lsiter + 1;
    else
        go=.false.
    endif
    
    if(my_rank.eq.0) then
	write(*,*)'      >lsiter, alpha, f, abs(fp-ft)/fp',lsiter, sngl(alpha), ft,abs(fp-ft)/fp
    endif
    if(ft.ge.fp) then
	if(abs(fp-ft)/fp.lt.f_tol) go=.false.
        alpha=(alpha+mu)/2;
    elseif(ft.lt.fp) then
	if(abs(fp-ft)/fp.lt.f_tol) go=.false.
        fp=ft;
        mu=alpha;
        alpha=2*alpha;
    endif
enddo
if(ft.gt.fp) then
	 alpha=mu;
endif
if(alpha.ne.0) then
	inv_x=x0+alpha*p;
	if(my_rank.eq.0) then
		! saving the updated velocity to disk
		do istg=1,size(inv_x)/nx_in/nz_in
			call filename(junk,trim(out_dir)//'/vel_band#',iband,'_stage#',istg,'_itr#',itr)
			call makesu2d(reshape(sngl(inv_x(nz_in*nx_in*(istg-1)+1:&
			       nz_in*nx_in*(istg))),(/nz_in,nx_in/))&
	                              ,junk,0.0,0.0,dx,dx,nx_in,nz_in,10)
		enddo	
	endif
	call grad_func(x=inv_x,gradient=inv_g,functional=inv_f)
else
	if(present(flag)) flag=.false.
endif
if(present(alpha_out)) alpha_out=alpha

end subroutine wWolfeLS
!********************************************************************************************

end module fwi_xeon


