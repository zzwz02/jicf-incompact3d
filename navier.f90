!################################################################################
!This file is part of Incompact3d.
!
!Incompact3d
!Copyright (c) 2012 Eric Lamballais and Sylvain Laizet
!eric.lamballais@univ-poitiers.fr / sylvain.laizet@gmail.com
!
!    Incompact3d is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation.
!
!    Incompact3d is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with the code.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!    We kindly request that you cite Incompact3d in your publications and 
!    presentations. The following citations are suggested:
!
!    1-Laizet S. & Lamballais E., 2009, High-order compact schemes for 
!    incompressible flows: a simple and efficient method with the quasi-spectral 
!    accuracy, J. Comp. Phys.,  vol 228 (15), pp 5989-6015
!
!    2-Laizet S. & Li N., 2011, Incompact3d: a powerful tool to tackle turbulence 
!    problems with up to 0(10^5) computational cores, Int. J. of Numerical 
!    Methods in Fluids, vol 67 (11), pp 1735-1757
!################################################################################

!********************************************************************
!
subroutine  intt (ux,uy,uz,gx,gy,gz,hx,hy,hz,ta1,tb1,tc1)
! 
!********************************************************************

USE param
USE variables
USE decomp_2d

implicit none

integer :: ijk,nxyz
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: gx,gy,gz
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: hx,hy,hz
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1

nxyz=xsize(1)*xsize(2)*xsize(3)

if ((nscheme.eq.1).or.(nscheme.eq.2)) then
if ((nscheme.eq.1.and.itime.eq.1.and.ilit.eq.0).or.&
     (nscheme.eq.2.and.itr.eq.1)) then
   do ijk=1,nxyz
      ux(ijk,1,1)=gdt(itr)*ta1(ijk,1,1)+ux(ijk,1,1)
      uy(ijk,1,1)=gdt(itr)*tb1(ijk,1,1)+uy(ijk,1,1) 
      uz(ijk,1,1)=gdt(itr)*tc1(ijk,1,1)+uz(ijk,1,1)
      gx(ijk,1,1)=ta1(ijk,1,1)
      gy(ijk,1,1)=tb1(ijk,1,1)
      gz(ijk,1,1)=tc1(ijk,1,1)            
   enddo
else
   if (nz.gt.1) then
      do ijk=1,nxyz
         ux(ijk,1,1)=adt(itr)*ta1(ijk,1,1)+bdt(itr)*gx(ijk,1,1)+ux(ijk,1,1)
         uy(ijk,1,1)=adt(itr)*tb1(ijk,1,1)+bdt(itr)*gy(ijk,1,1)+uy(ijk,1,1)   
         uz(ijk,1,1)=adt(itr)*tc1(ijk,1,1)+bdt(itr)*gz(ijk,1,1)+uz(ijk,1,1)
         gx(ijk,1,1)=ta1(ijk,1,1)
         gy(ijk,1,1)=tb1(ijk,1,1)
         gz(ijk,1,1)=tc1(ijk,1,1)            
      enddo
   else
      do ijk=1,nxyz
         ux(ijk,1,1)=adt(itr)*ta1(ijk,1,1)+bdt(itr)*gx(ijk,1,1)+ux(ijk,1,1)
         uy(ijk,1,1)=adt(itr)*tb1(ijk,1,1)+bdt(itr)*gy(ijk,1,1)+uy(ijk,1,1)   
         gx(ijk,1,1)=ta1(ijk,1,1)
         gy(ijk,1,1)=tb1(ijk,1,1)
      enddo
   endif
endif
endif

if (nscheme.eq.3) then 
   if (nz.gt.1) then
      if (adt(itr)==0.) then
         do ijk=1,nxyz
            gx(ijk,1,1)=dt*ta1(ijk,1,1)
            gy(ijk,1,1)=dt*tb1(ijk,1,1)
            gz(ijk,1,1)=dt*tc1(ijk,1,1)
         enddo
      else
         do ijk=1,nxyz
            gx(ijk,1,1)=adt(itr)*gx(ijk,1,1)+dt*ta1(ijk,1,1)
            gy(ijk,1,1)=adt(itr)*gy(ijk,1,1)+dt*tb1(ijk,1,1)
            gz(ijk,1,1)=adt(itr)*gz(ijk,1,1)+dt*tc1(ijk,1,1)
         enddo
      endif
      do ijk=1,nxyz
         ux(ijk,1,1)=ux(ijk,1,1)+bdt(itr)*gx(ijk,1,1)
         uy(ijk,1,1)=uy(ijk,1,1)+bdt(itr)*gy(ijk,1,1)
         uz(ijk,1,1)=uz(ijk,1,1)+bdt(itr)*gz(ijk,1,1)
      enddo
   else
      if (adt(itr)==0.) then
         do ijk=1,nxyz
            gx(ijk,1,1)=dt*ta1(ijk,1,1)
            gy(ijk,1,1)=dt*tb1(ijk,1,1)
         enddo
      else
         do ijk=1,nxyz
            gx(ijk,1,1)=adt(itr)*gx(ijk,1,1)+dt*ta1(ijk,1,1)
            gy(ijk,1,1)=adt(itr)*gy(ijk,1,1)+dt*tb1(ijk,1,1)
         enddo
      endif
      do ijk=1,nxyz
         ux(ijk,1,1)=ux(ijk,1,1)+bdt(itr)*gx(ijk,1,1)
         uy(ijk,1,1)=uy(ijk,1,1)+bdt(itr)*gy(ijk,1,1)
      enddo
   endif
endif

if (nscheme==4) then
   if ((itime.eq.1).and.(ilit.eq.0)) then
      if (nrank==0) print *,'start with Euler',itime
      do ijk=1,nxyz !start with Euler
         ux(ijk,1,1)=dt*ta1(ijk,1,1)+ux(ijk,1,1)
         uy(ijk,1,1)=dt*tb1(ijk,1,1)+uy(ijk,1,1) 
         uz(ijk,1,1)=dt*tc1(ijk,1,1)+uz(ijk,1,1)
         gx(ijk,1,1)=ta1(ijk,1,1)
         gy(ijk,1,1)=tb1(ijk,1,1)
         gz(ijk,1,1)=tc1(ijk,1,1)            
      enddo
   else
      if  ((itime.eq.2).and.(ilit.eq.0)) then
          if (nrank==0) print *,'then with AB2',itime
         do ijk=1,nxyz
            ux(ijk,1,1)=1.5*dt*ta1(ijk,1,1)-0.5*dt*gx(ijk,1,1)+ux(ijk,1,1)
            uy(ijk,1,1)=1.5*dt*tb1(ijk,1,1)-0.5*dt*gy(ijk,1,1)+uy(ijk,1,1)
            uz(ijk,1,1)=1.5*dt*tc1(ijk,1,1)-0.5*dt*gz(ijk,1,1)+uz(ijk,1,1)
            hx(ijk,1,1)=gx(ijk,1,1)
            hy(ijk,1,1)=gy(ijk,1,1)
            hz(ijk,1,1)=gz(ijk,1,1)
            gx(ijk,1,1)=ta1(ijk,1,1)
            gy(ijk,1,1)=tb1(ijk,1,1)
            gz(ijk,1,1)=tc1(ijk,1,1)
         enddo 
      else
         do ijk=1,nxyz
            ux(ijk,1,1)=adt(itr)*ta1(ijk,1,1)+bdt(itr)*gx(ijk,1,1)+&
                 cdt(itr)*hx(ijk,1,1)+ux(ijk,1,1)
            uy(ijk,1,1)=adt(itr)*tb1(ijk,1,1)+bdt(itr)*gy(ijk,1,1)+&
                 cdt(itr)*hy(ijk,1,1)+uy(ijk,1,1)
            uz(ijk,1,1)=adt(itr)*tc1(ijk,1,1)+bdt(itr)*gz(ijk,1,1)+&
                 cdt(itr)*hz(ijk,1,1)+uz(ijk,1,1)
            hx(ijk,1,1)=gx(ijk,1,1)
            hy(ijk,1,1)=gy(ijk,1,1)
            hz(ijk,1,1)=gz(ijk,1,1)
            gx(ijk,1,1)=ta1(ijk,1,1)
            gy(ijk,1,1)=tb1(ijk,1,1)
            gz(ijk,1,1)=tc1(ijk,1,1)
         enddo
      endif
   endif
endif


return
end subroutine intt

!********************************************************************
!
subroutine corgp (ux,gx,uy,uz,px,py,pz,phi)
! 
!********************************************************************

USE decomp_2d
USE variables
USE param
USE var
USE MPI

implicit none

integer :: ijk,nxyz
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,px,py,pz,phi
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: gx,phi22

nxyz=xsize(1)*xsize(2)*xsize(3)

do ijk=1,nxyz
   ux(ijk,1,1)=-px(ijk,1,1)+ux(ijk,1,1)
   uy(ijk,1,1)=-py(ijk,1,1)+uy(ijk,1,1)
   uz(ijk,1,1)=-pz(ijk,1,1)+uz(ijk,1,1)
enddo

if (itype==2) then !channel flow
   call transpose_x_to_y(ux,gx)
   if (iscalar==1) call transpose_x_to_y(phi,phi22)
   call channel(gx,phi22)
   if (iscalar==1) call transpose_y_to_x(phi22,phi)
   call transpose_y_to_x(gx,ux)
endif

return
end subroutine corgp

!*********************************************************
!
subroutine inflow (ux,uy,uz,phi)
!  
!*********************************************************

USE param
USE IBM
USE variables
USE decomp_2d

implicit none

integer  :: k,j
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,phi
real(mytype) :: r1,r2,r3,y,um



if (ivirt==0) then
  ibm_y_max=1
else
   ibm_y_max=97
  !ibm_y_max=min((itime-0)/5000+91,97)
endif

call ecoule(ux,uy,uz)
call random_number(bxo)
call random_number(byo)
call random_number(bzo)

if (itype.ne.5) then
   if (iin.eq.1) then  
      do k=1,xsize(3)
      do j=1,xsize(2)
         bxx1(j,k)=bxx1(j,k)+bxo(j,k)*noise1
         bxy1(j,k)=bxy1(j,k)+byo(j,k)*noise1
         bxz1(j,k)=bxz1(j,k)+bzo(j,k)*noise1
         !if (print_flag==1 .and. (k+xstart(3)-1)==1) print *,'bxx1',(j+xstart(2)-1),k,bxx1(j,k)
      enddo
      enddo
      if (iscalar==1) then
         do k=1,xsize(3)
         do j=1,xsize(2)
            phi(1,j,k)=0.
         enddo
         enddo
      endif
   endif
elseif (itype==5) then
   if (rezone>0) then
      if (iin.eq.1) then  
         do k=1,xsize(3)
         do j=1,xsize(2)
            bxx1(j,k)=ux(rezone,j,k)
            bxy1(j,k)=uy(rezone,j,k)
            bxz1(j,k)=uz(rezone,j,k)
         enddo
         enddo
         if (iscalar==1) then
            do k=1,xsize(3)
            do j=1,xsize(2)
               !phi(1,j,k)=phi(rezone,j,k)
               phi(1,j,k)=0.
            enddo
            enddo
         endif !(iscalar==1)
      endif !(iin.eq.1)
   else  !rezone==0, no recycling

      if (iin.eq.1) then  
         do k=1,xsize(3)
         do j=1,xsize(2)
	    bxx1(j,k)=bxx1(j,k)+bxo(j,k)*noise1
      	    bxy1(j,k)=bxy1(j,k)+byo(j,k)*noise1
      	    bxz1(j,k)=bxz1(j,k)+bzo(j,k)*noise1
            !if (print_flag==1 .and. (k+xstart(3)-1)==1) print *,'bxx1',(j+xstart(2)-1),k,bxx1(j,k)
         enddo
         enddo
         if (iscalar==1) then
            do k=1,xsize(3)
            do j=1,xsize(2)
               phi(1,j,k)=0.
            enddo
            enddo
         endif !(iscalar==1)
      endif !(iin.eq.1)
   endif !(rezone==0)
endif !(itype.==.5)

return
end subroutine inflow 

!*********************************************************
!
subroutine outflow (ux,uy,uz,phi)
!
!*********************************************************

USE param
USE variables
USE decomp_2d
USE MPI

implicit none

integer :: j,k,i, code
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,phi
real(mytype) :: udx,udy,udz,uddx,uddy,uddz,uxmax,&
     uxmin,vphase,cx,coef,uxmax1,uxmin1



udx=1./dx
udy=1./dy
udz=1./dz
uddx=0.5/dx
uddy=0.5/dy
uddz=0.5/dz
cx=0.5*(u1+u2)*gdt(itr)*udx


uxmax=-1609.
uxmin=1609.
do k=1,xsize(3)
do j=1,xsize(2)
   if (ux(nx-1,j,k).gt.uxmax) uxmax=ux(nx-1,j,k)
   if (ux(nx-1,j,k).lt.uxmin) uxmin=ux(nx-1,j,k)
enddo
enddo
call MPI_ALLREDUCE(uxmax,uxmax1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
call MPI_ALLREDUCE(uxmin,uxmin1,1,real_type,MPI_MIN,MPI_COMM_WORLD,code)
vphase=0.5*(uxmax1+uxmin1)
cx=vphase*gdt(itr)*udx


if (itype.ne.9) then
   do k=1,xsize(3)
   do j=1,xsize(2)
    if ((j+xstart(2)-1)<ibm_y_max) then
	bxxn(j,k)=0.
      	bxyn(j,k)=0.
      	bxzn(j,k)=0.
     else
     	bxxn(j,k)=ux(nx,j,k)-cx*(ux(nx,j,k)-ux(nx-1,j,k))
      	bxyn(j,k)=uy(nx,j,k)-cx*(uy(nx,j,k)-uy(nx-1,j,k))
      	bxzn(j,k)=uz(nx,j,k)-cx*(uz(nx,j,k)-uz(nx-1,j,k))
     endif
      !bxxn(j,k)=ux(nx,j,k)-cx*(ux(nx,j,k)-ux(nx-1,j,k))
      !bxyn(j,k)=uy(nx,j,k)-cx*(uy(nx,j,k)-uy(nx-1,j,k))
      !bxzn(j,k)=uz(nx,j,k)-cx*(uz(nx,j,k)-uz(nx-1,j,k))
   enddo
   enddo
   if (iscalar==1) then
      do k=1,xsize(3)
      do j=1,xsize(2)
         phi(nx,j,k)=phi(nx,j,k)-cx*(phi(nx,j,k)-phi(nx-1,j,k))
      enddo
      enddo
   endif
else
print *,'NOT READY'
stop
endif


return
end subroutine outflow 

!**********************************************************************
!
subroutine ecoule(ux1,uy1,uz1)
!
!**********************************************************************

USE param
USE IBM
USE variables
USE decomp_2d

implicit none

integer  :: i,j,k,jj1,jj2 
real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
real(mytype) :: x,y,z,ym,a,b
real(mytype) :: r1,r2,r3,r
real(mytype) :: uh,ud,um,xv,bruit1

bxx1=0.;bxy1=0.;bxz1=0.
byx1=0.;byy1=0.;byz1=0.
bzx1=0.;bzy1=0.;bzz1=0. 

!ITYPE=1 --> Constant flow field
!ITYPE=2 --> Channel flow
!ITYPE=3 --> Wake flow
!ITYPE=4 --> Mixing layer with splitter plate
!ITYPE=5 --> Jet in cross flow
!ITYPE=6 --> Taylor Green vortices
!ITYPE=7 --> Cavity flow
!ITYPE=8 --> Flat plate Boundary layer
!ITYPE=9 --> Tank 

if (itype.eq.1) then
   um=0.5*(u1+u2)
   do k=1,xsize(3)
   do j=1,xsize(2)
      bxx1(j,k)=um
      bxy1(j,k)=0.
      bxz1(j,k)=0.
   enddo
   enddo
endif

if (itype.eq.2) then
   do k=1,xsize(3)
   do j=1,xsize(2)
      if (istret.eq.0) y=(j+xstart(2)-1-1)*dy
      if (istret.ne.0) y=yp(j+xstart(2)-1)
!      print *,nrank,j+xstart(2)-1,yp(j+xstart(2)-1),1.-y*y
         bxx1(j,k)=-1.5*y*(y-2)  !-6./(a-b)^2*(y-a)*(y-b)*v_m
   enddo
   
   enddo 
endif

if (itype.eq.3) then
   
endif

if (itype.eq.4) then
  
endif

if (itype.eq.5) then
   if (nclx.eq.0) then
      print *,'X NOT POSSIBLE'
      stop
   endif
   if (ncly==0) then
      print *,'Y NOT POSSIBLE'
      stop
   endif
   if (nclz.ne.0) then
      print *,'Z NOT POSSIBLE'
      stop
   endif

   do k=1,xsize(3)
   do j=1,xsize(2)
     if (j+xstart(2)-1<ibm_y_max) then
	bxx1(j,k)=0
     else
      if (istret.eq.0) y=(j+xstart(2)-1-1)*dy
      if (istret.ne.0) y=yp(j+xstart(2)-1)
!      print *,nrank,j+xstart(2)-1,yp(j+xstart(2)-1),1.-y*y
	a=yly
	b=yp(ibm_y_max)
        bxx1(j,k)=-6*(y-a)*(y-b)/((a-b)*(a-b))
      !bxx1(j,k)=-1.5*(y-2.)*(y-4.)
     endif

   enddo
   enddo
endif

if (itype.eq.6) then
   t=0.
   xv=1./100.
   xxk1=twopi/xlx
   xxk2=twopi/yly
   do k=1,xsize(3)
      z=(k+xstart(3)-1-1)*dz
   do j=1,xsize(2)
      y=(j+xstart(2)-1-1)*dy
      do i=1,xsize(1)
         x=(i-1)*dx
         ux1(i,j,k)=sin(2.*pi*x)*cos(2.*pi*y)*cos(2.*pi*z)
         uy1(i,j,k)=sin(2.*pi*y)*cos(2.*pi*x)*cos(2.*pi*z)
         uz1(i,j,k)=sin(2.*pi*z)*cos(2.*pi*x)*cos(2.*pi*y)
         bxx1(j,k)=0.
         bxy1(j,k)=0.
         bxz1(j,k)=0.
      enddo
   enddo
   enddo   
endif

if (itype.eq.7) then

endif

if (itype.eq.8) then
   
endif

if (itype.eq.9) then
  
endif


if (itype.eq.10) then
   do k=1,xsize(3)
   do j=1,xsize(2)
      bxx1(j,k)=0.
      bxy1(j,k)=0.
      bxz1(j,k)=0.
   enddo
   enddo  
endif

return
end subroutine ecoule

!********************************************************************
!
subroutine init (ux1,uy1,uz1,ep1,phi1,gx1,gy1,gz1,phis1,hx1,hy1,hz1,phiss1)
!
!********************************************************************

USE decomp_2d
USE decomp_2d_io
USE variables
USE param
USE MPI

implicit none

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,phi1,ep1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: gx1,gy1,gz1,phis1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: hx1,hy1,hz1,phiss1

real(mytype) :: y,r,um,r1,r2,r3
integer :: k,j,i,fh,ierror,code,ii
integer (kind=MPI_OFFSET_KIND) :: disp

if (iin.eq.1) then !generation of a random noise


   call system_clock(count=code)
   call random_seed(size = ii)
   call random_seed(put = code+63949*nrank*(/ (i - 1, i = 1, ii) /)) ! 63949, "large" prime number

   call random_number(ux1)
   call random_number(uy1)
   call random_number(uz1)



   do k=1,xsize(3)
   do j=1,xsize(2)
   do i=1,xsize(1)
      ux1(i,j,k)=noise*(ux1(i,j,k)-0.5)
      uy1(i,j,k)=noise*(uy1(i,j,k)-0.5)
      uz1(i,j,k)=noise*(uz1(i,j,k)-0.5)
   enddo
   enddo
   enddo

!modulation of the random noise
   do k=1,xsize(3)
   do j=1,xsize(2)
      if (istret.eq.0) y=(j+xstart(2)-1-1)*dy-yly/2.
      if (istret.ne.0) y=yp(j+xstart(2)-1)-yly/2.
      um=exp(-0.2*y*y)
      do i=1,xsize(1)
         ux1(i,j,k)=um*ux1(i,j,k)
         uy1(i,j,k)=um*uy1(i,j,k)
         uz1(i,j,k)=um*uz1(i,j,k)
      enddo
   enddo
   enddo

   if (iscalar==1) then
      do k=1,xsize(3)
      do j=1,xsize(2)
      do i=1,xsize(1)
   !      if ((j+xstart(2)-1).ge.nym) then
   !         phi1(i,j,k)=1.
   !      else
         phi1(i,j,k)=0.
   !      endif
         phis1(i,j,k)=phi1(i,j,k)
         phiss1(i,j,k)=phis1(i,j,k)
      enddo
      enddo
      enddo
   endif
endif

if (iin.eq.2) then !read a correlated noise generated before
endif

!MEAN FLOW PROFILE
call ecoule(ux1,uy1,uz1)
!INIT FOR G AND U=MEAN FLOW + NOISE
do k=1,xsize(3)
do j=1,xsize(2)
do i=1,xsize(1)
   ux1(i,j,k)=ux1(i,j,k)+bxx1(j,k)
   uy1(i,j,k)=uy1(i,j,k)+bxy1(j,k)
   uz1(i,j,k)=uz1(i,j,k)+bxz1(j,k)
   gx1(i,j,k)=ux1(i,j,k)
   gy1(i,j,k)=uy1(i,j,k)
   gz1(i,j,k)=uz1(i,j,k)
   hx1(i,j,k)=gx1(i,j,k)
   hy1(i,j,k)=gy1(i,j,k)
   hz1(i,j,k)=gz1(i,j,k)
enddo
enddo
enddo

ibm_y_max=1
if (ivirt==2) then
   call MPI_FILE_OPEN(MPI_COMM_WORLD, 'epsilon.dat', &
        MPI_MODE_RDONLY, MPI_INFO_NULL, &
        fh, ierror)
   disp = 0_MPI_OFFSET_KIND
   call decomp_2d_read_var(fh,disp,1,ep1) 
   call MPI_FILE_CLOSE(fh,ierror)
   if (nrank==0) print *,'read epsilon file done from init'
   print *,ep1
endif


return
end subroutine init

!********************************************************************
!
subroutine divergence (ux1,uy1,uz1,ep1,ta1,tb1,tc1,di1,td1,te1,tf1,&
     td2,te2,tf2,di2,ta2,tb2,tc2,ta3,tb3,tc3,di3,td3,te3,tf3,pp3,&
     nxmsize,nymsize,nzmsize,ph1,ph3,ph4,nlock)
! 
!********************************************************************

USE param
USE IBM
USE decomp_2d
USE variables
USE MPI
   
implicit none

TYPE(DECOMP_INFO) :: ph1,ph3,ph4

integer :: nxmsize,nymsize,nzmsize

!X PENCILS NX NY NZ  -->NXM NY NZ
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,di1,ux1,uy1,uz1,ep1
real(mytype),dimension(nxmsize,xsize(2),xsize(3)) :: td1,te1,tf1 
!Y PENCILS NXM NY NZ  -->NXM NYM NZ
real(mytype),dimension(ph1%yst(1):ph1%yen(1),ysize(2),ysize(3)) :: td2,te2,tf2,di2
real(mytype),dimension(ph1%yst(1):ph1%yen(1),nymsize,ysize(3)) :: ta2,tb2,tc2
!Z PENCILS NXM NYM NZ  -->NXM NYM NZM
real(mytype),dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),zsize(3)) :: ta3,tb3,tc3,di3
real(mytype),dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize) :: td3,te3,tf3,pp3



integer :: ijk,nvect1,nvect2,nvect3,i,j,k,nlock
integer :: code
real(mytype) :: tmax,tmoy,tmax1,tmoy1


nvect1=xsize(1)*xsize(2)*xsize(3)
nvect2=ysize(1)*ysize(2)*ysize(3)
nvect3=(ph1%zen(1)-ph1%zst(1)+1)*(ph1%zen(2)-ph1%zst(2)+1)*nzmsize

if (nlock==1) then
   if (ivirt.eq.0) ep1(:,:,:)=0.
   do ijk=1,nvect1
      ta1(ijk,1,1)=(1.-ep1(ijk,1,1))*ux1(ijk,1,1)
      tb1(ijk,1,1)=(1.-ep1(ijk,1,1))*uy1(ijk,1,1)
      tc1(ijk,1,1)=(1.-ep1(ijk,1,1))*uz1(ijk,1,1)
   enddo
else
   ta1(:,:,:)=ux1(:,:,:)
   tb1(:,:,:)=uy1(:,:,:)
   tc1(:,:,:)=uz1(:,:,:)
endif



!WORK X-PENCILS
call decx6(td1,ta1,di1,sx,cfx6,csx6,cwx6,xsize(1),nxmsize,xsize(2),xsize(3),0)
call inter6(te1,tb1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)
call inter6(tf1,tc1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)

call transpose_x_to_y(td1,td2,ph4)!->NXM NY NZ
call transpose_x_to_y(te1,te2,ph4)
call transpose_x_to_y(tf1,tf2,ph4)


!WORK Y-PENCILS
call intery6(ta2,td2,di2,sy,cifyp6,cisyp6,ciwyp6,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),1)
call decy6(tb2,te2,di2,sy,cfy6,csy6,cwy6,ppyi,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),0)
call intery6(tc2,tf2,di2,sy,cifyp6,cisyp6,ciwyp6,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),1)

call transpose_y_to_z(ta2,ta3,ph3)!->NXM NYM NZ
call transpose_y_to_z(tb2,tb3,ph3)
call transpose_y_to_z(tc2,tc3,ph3)


!WORK Z-PENCILS
call interz6(td3,ta3,di3,sz,cifzp6,ciszp6,ciwzp6,(ph1%zen(1)-ph1%zst(1)+1),&
     (ph1%zen(2)-ph1%zst(2)+1),zsize(3),nzmsize,1)    
call interz6(te3,tb3,di3,sz,cifzp6,ciszp6,ciwzp6,(ph1%zen(1)-ph1%zst(1)+1),&
     (ph1%zen(2)-ph1%zst(2)+1),zsize(3),nzmsize,1)
call decz6(tf3,tc3,di3,sz,cfz6,csz6,cwz6,(ph1%zen(1)-ph1%zst(1)+1),&
     (ph1%zen(2)-ph1%zst(2)+1),zsize(3),nzmsize,0)


do k=1,nzmsize
do j=ph1%zst(2),ph1%zen(2)
do i=ph1%zst(1),ph1%zen(1)
pp3(i,j,k)=td3(i,j,k)+te3(i,j,k)+tf3(i,j,k)
enddo
enddo
enddo


if (nlock==2) then
   pp3(:,:,:)=pp3(:,:,:)-pp3(ph1%zst(1),ph1%zst(2),nzmsize)
endif

tmax=-1609.
tmoy=0.
do k=1,nzmsize
do j=ph1%zst(2),ph1%zen(2)
do i=ph1%zst(1),ph1%zen(1)
   if (pp3(i,j,k).gt.tmax) tmax=pp3(i,j,k)     
   tmoy=tmoy+abs(pp3(i,j,k))
enddo
enddo
enddo
tmoy=tmoy/nvect3

call MPI_REDUCE(tmax,tmax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
call MPI_REDUCE(tmoy,tmoy1,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)!

if (print_flag==1 .and. nrank==0) then
     if (nlock==2) then
        print *,'DIV U final Max=',tmax1
        print *,'DIV U final Moy=',tmoy1/real(nproc)
     else
        print *,'DIV U* Max=',tmax1
        print *,'DIV U* Moyy=',tmoy1/real(nproc)
     endif
endif

return
end subroutine divergence

!********************************************************************
!
subroutine gradp(ta1,tb1,tc1,di1,td2,tf2,ta2,tb2,tc2,di2,&
     ta3,tc3,di3,pp3,nxmsize,nymsize,nzmsize,ph2,ph3)
!
!********************************************************************

USE param 
USE decomp_2d
USE variables

implicit none

TYPE(DECOMP_INFO) :: ph2,ph3
integer :: i,j,k,ijk,nxmsize,nymsize,nzmsize,code
integer, dimension(2) :: dims, dummy_coords
logical, dimension(2) :: dummy_periods

real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),nzmsize) :: pp3 
!Z PENCILS NXM NYM NZM-->NXM NYM NZ
real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),zsize(3)) :: ta3,tc3,di3 
!Y PENCILS NXM NYM NZ -->NXM NY NZ
real(mytype),dimension(ph3%yst(1):ph3%yen(1),nymsize,ysize(3)) :: ta2,tc2
real(mytype),dimension(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)) :: tb2,td2,tf2,di2
!X PENCILS NXM NY NZ  -->NX NY NZ
real(mytype),dimension(nxmsize,xsize(2),xsize(3)) :: td1,te1,tf1 
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,di1 




!WORK Z-PENCILS

call interiz6(ta3,pp3,di3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
     (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
call deciz6(tc3,pp3,di3,sz,cfip6z,csip6z,cwip6z,cfz6,csz6,cwz6,&
     (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)

!WORK Y-PENCILS
call transpose_z_to_y(ta3,ta2,ph3) !nxm nym nz
call transpose_z_to_y(tc3,tc2,ph3)

call interiy6(tb2,ta2,di2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
     (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
call deciy6(td2,ta2,di2,sy,cfip6y,csip6y,cwip6y,cfy6,csy6,cwy6,ppy,&
     (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
call interiy6(tf2,tc2,di2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
     (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)

!WORK X-PENCILS

call transpose_y_to_x(tb2,td1,ph2) !nxm ny nz
call transpose_y_to_x(td2,te1,ph2)
call transpose_y_to_x(tf2,tf1,ph2)

call deci6(ta1,td1,di1,sx,cfip6,csip6,cwip6,cfx6,csx6,cwx6,&
     nxmsize,xsize(1),xsize(2),xsize(3),1)
call interi6(tb1,te1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
     nxmsize,xsize(1),xsize(2),xsize(3),1)
call interi6(tc1,tf1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
     nxmsize,xsize(1),xsize(2),xsize(3),1)


!we are in X pencils:
do k=1,xsize(3)
do j=1,xsize(2)
   dpdyx1(j,k)=tb1(1,j,k)/gdt(itr)
   dpdzx1(j,k)=tc1(1,j,k)/gdt(itr)
   dpdyxn(j,k)=tb1(xsize(1),j,k)/gdt(itr)
   dpdzxn(j,k)=tc1(xsize(1),j,k)/gdt(itr)
enddo
enddo


if (xstart(3)==1) then
   do j=1,xsize(2)
   do i=1,xsize(1)
      dpdxz1(i,j)=ta1(i,j,1)/gdt(itr)
      dpdyz1(i,j)=tb1(i,j,1)/gdt(itr)
   enddo
   enddo
endif
if (xend(3)==nz) then
   do j=1,xsize(2)
   do i=1,xsize(1)
      dpdxzn(i,j)=ta1(i,j,xsize(3))/gdt(itr)
      dpdyzn(i,j)=tb1(i,j,xsize(3))/gdt(itr)
   enddo
   enddo
endif


   ! determine the processor grid in use
   call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, &
         dims, dummy_periods, dummy_coords, code)

   if (dims(1)==1) then
      do k=1,xsize(3)
      do i=1,xsize(1)
      dpdxy1(i,k)=ta1(i,1,k)/gdt(itr)
      dpdzy1(i,k)=tc1(i,1,k)/gdt(itr)
      enddo
      enddo
      do k=1,xsize(3)
      do i=1,xsize(1)
      dpdxyn(i,k)=ta1(i,xsize(2),k)/gdt(itr)
      dpdzyn(i,k)=tc1(i,xsize(2),k)/gdt(itr)
      enddo
      enddo
   else
!find j=1 and j=ny
      if (xstart(2)==1) then
         do k=1,xsize(3)
         do i=1,xsize(1)
      dpdxy1(i,k)=ta1(i,1,k)/gdt(itr)
      dpdzy1(i,k)=tc1(i,1,k)/gdt(itr)
         enddo
         enddo
      endif
!      print *,nrank,xstart(2),ny-(nym/p_row)
       if (ny-(nym/dims(1))==xstart(2)) then
         do k=1,xsize(3)
         do i=1,xsize(1)
      dpdxyn(i,k)=ta1(i,xsize(2),k)/gdt(itr)
      dpdzyn(i,k)=tc1(i,xsize(2),k)/gdt(itr)
         enddo
         enddo
      endif

   endif

return
end subroutine gradp

!********************************************************************
!
subroutine corgp_IBM (ux,uy,uz,px,py,pz,nlock)
! 
!********************************************************************

USE param 
USE decomp_2d
USE variables

implicit none

integer :: ijk,nlock,nxyz
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,px,py,pz

if (itime==1) then
   px(:,:,:)=0.
   py(:,:,:)=0.
   pz(:,:,:)=0.
endif

nxyz=xsize(1)*xsize(2)*xsize(3)

if (nlock.eq.1) then
   if (nz.gt.1) then
      do ijk=1,nxyz
         uy(ijk,1,1)=-py(ijk,1,1)+uy(ijk,1,1) 
         uz(ijk,1,1)=-pz(ijk,1,1)+uz(ijk,1,1) 
         ux(ijk,1,1)=-px(ijk,1,1)+ux(ijk,1,1)
      enddo
   else
      do ijk=1,nxyz
         uy(ijk,1,1)=-py(ijk,1,1)+uy(ijk,1,1) 
         ux(ijk,1,1)=-px(ijk,1,1)+ux(ijk,1,1)
      enddo
   endif
endif
if (nlock.eq.2) then
   if (nz.gt.1) then
      do ijk=1,nxyz
         uy(ijk,1,1)=py(ijk,1,1)+uy(ijk,1,1) 
         uz(ijk,1,1)=pz(ijk,1,1)+uz(ijk,1,1) 
         ux(ijk,1,1)=px(ijk,1,1)+ux(ijk,1,1)
      enddo
   else
      do ijk=1,nxyz
         uy(ijk,1,1)=py(ijk,1,1)+uy(ijk,1,1) 
         ux(ijk,1,1)=px(ijk,1,1)+ux(ijk,1,1)
      enddo
   endif
endif

return
end subroutine corgp_IBM

!*******************************************************************
!
subroutine body(ux,uy,uz,ep1)
!
!*******************************************************************

USE param 
USE decomp_2d
USE variables
USE IBM
!USE decomp_2d_io

implicit none

real(mytype), dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ux,uy,uz,ep1
integer :: i,j,k,tempa,tempb
real(mytype) :: r2
!real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu
!character(len=20) :: filename

if ((nrank==0) .and. (print_flag==1)) print *,"ibm_y_max",ibm_y_max

!if (itime==ifirst) then
if (.true.) then
  ep1=0.
  do k=xstart(3),xend(3)
  do j=xstart(2),xend(2)
  do i=xstart(1),xend(1)
    tempa=xjicf-i
    tempb=nz/2-k
    r2=dx2*real(tempa*tempa)+dz2*real(tempb*tempb)
  
    if ((j .le. ibm_y_max) .and. (r2 > rjicf)) then
      ux(i,j,k)=0.
      uy(i,j,k)=0.
      uz(i,j,k)=0.
      ep1(i,j,k)=1.
    endif
  enddo
  enddo
  enddo
else
  do k=xstart(3),xend(3)
  do j=xstart(2),xend(2)
  do i=xstart(1),xend(1)
    if (ep1(i,j,k)==1.) then
      ux(i,j,k)=0.
      uy(i,j,k)=0.
      uz(i,j,k)=0.
    endif
  enddo
  enddo
  enddo
endif

!if (mod(itime,imodulo)==0) then
!  uvisu=0.
!  call fine_to_coarseV(1,ep1,uvisu)
!  995 format('ibm',I3.3)
!      write(filename, 995) itime/imodulo
!  call decomp_2d_write_one(1,uvisu,filename,2)
!endif

return  
end subroutine body

!****************************************************************************
!
subroutine pre_correc(ux,uy,uz)
!
!****************************************************************************

USE decomp_2d
USE variables
USE param
USE var
USE MPI


implicit none

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: gx
integer :: i,j,k,code,njicf1,tempa,tempb
real(mytype) :: ut,ut1,utt,ut11,ujicf1,r2,temp1,temp
integer, dimension(2) :: dims, dummy_coords
logical, dimension(2) :: dummy_periods

if (itime==1) then
   dpdyx1=0.
   dpdzx1=0.
   dpdyxn=0.
   dpdzxn=0.
   dpdxy1=0.
   dpdzy1=0.
   dpdxyn=0.
   dpdzyn=0.
   dpdxz1=0.
   dpdyz1=0.
   dpdxzn=0.
   dpdyzn=0.
endif


!we are in X pencils:
do k=1,xsize(3)
do j=1,xsize(2)
   dpdyx1(j,k)=dpdyx1(j,k)*gdt(itr)
   dpdzx1(j,k)=dpdzx1(j,k)*gdt(itr)
   dpdyxn(j,k)=dpdyxn(j,k)*gdt(itr)
   dpdzxn(j,k)=dpdzxn(j,k)*gdt(itr)
enddo
enddo

if (xstart(3)==1) then
   do j=1,xsize(2)
   do i=1,xsize(1)
      dpdxz1(i,j)=dpdxz1(i,j)*gdt(itr)
      dpdyz1(i,j)=dpdyz1(i,j)*gdt(itr)
   enddo
   enddo
endif
if (xend(3)==nz) then
   do j=1,xsize(2)
   do i=1,xsize(1)
      dpdxzn(i,j)=dpdxzn(i,j)*gdt(itr)
      dpdyzn(i,j)=dpdyzn(i,j)*gdt(itr)
   enddo
   enddo
endif

if (xstart(2)==1) then
   do k=1,xsize(3)
   do i=1,xsize(1)
      dpdxy1(i,k)=dpdxy1(i,k)*gdt(itr)
      dpdzy1(i,k)=dpdzy1(i,k)*gdt(itr)
   enddo
   enddo
endif
if (xend(2)==ny) then
   do k=1,xsize(3)
   do i=1,xsize(1)
      dpdxyn(i,k)=dpdxyn(i,k)*gdt(itr)
      dpdzyn(i,k)=dpdzyn(i,k)*gdt(itr)
   enddo
   enddo
endif

!Computatation of the flow rate Inflow/Outflow
!we are in X pencils:
if (nclx==2) then
if (itype .eq. 2) then
   ut1=0.
   do k=1,xsize(3)
   do j=1,xsize(2)
      ut1=ut1+bxx1(j,k)
   enddo
   enddo
   ut1=ut1/xsize(2)/xsize(3)
   call MPI_ALLREDUCE(ut1,ut11,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
   ut11=ut11/nproc
   ut=0.
   do k=1,xsize(3)
   do j=1,xsize(2)
      ut=ut+bxxn(j,k)
   enddo
   enddo
   ut=ut/xsize(2)/xsize(3)
   call MPI_ALLREDUCE(ut,utt,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
   utt=utt/nproc

   if (print_flag==1 .and. nrank==0) print *,'FLOW RATE I/O',ut11,utt

   do k=1,xsize(3)
   do j=1,xsize(2)
      bxxn(j,k)=bxxn(j,k)-utt+ut11
   enddo
   enddo
else if (itype .eq. 5) then

   bxx11(:,:,:)=0
   bxxn1(:,:,:)=0
   do k=1,xsize(3)
   do j=1,xsize(2)
      bxx11(1,j,k)=bxx1(j,k)
      bxxn1(1,j,k)=bxxn(j,k)
   enddo
   enddo

   call transpose_x_to_y(ux,gx)
   call transpose_x_to_y(bxx11,bxx12)
   call transpose_x_to_y(bxxn1,bxxn2)

   if (rezone>0) call recycling(bxx12,gx)
   call jicf(bxx12,bxxn2)

   call transpose_y_to_x(gx,ux)
   call transpose_y_to_x(bxxn2,bxxn1)
   call transpose_y_to_x(bxx12,bxx11)

   do k=1,xsize(3)
   do j=1,xsize(2)
      bxx1(j,k)=bxx11(1,j,k)
      bxxn(j,k)=bxxn1(1,j,k)
   enddo
   enddo

endif
endif

!********NCLX==2*************************************
!****************************************************
if (nclx.eq.2) then
   do k=1,xsize(3)
   do j=1,xsize(2)
      ux(1 ,j,k)=bxx1(j,k)
      uy(1 ,j,k)=bxy1(j,k)+dpdyx1(j,k)
      uz(1 ,j,k)=bxz1(j,k)+dpdzx1(j,k)
      ux(nx,j,k)=bxxn(j,k)
      uy(nx,j,k)=bxyn(j,k)+dpdyxn(j,k)
      uz(nx,j,k)=bxzn(j,k)+dpdzxn(j,k)
   enddo
   enddo
endif
!****************************************************
!********NCLY==2*************************************
!****************************************************
!WE ARE IN X PENCIL!!!!!!
if ((ncly==2) .and. (iimplicit==0))then
   if (itype.eq.2) then

   ! determine the processor grid in use
   call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, &
         dims, dummy_periods, dummy_coords, code)

   if (dims(1)==1) then
      do k=1,xsize(3)
      do i=1,xsize(1)
         ux(i,1,k)=0.+dpdxy1(i,k)
         uy(i,1,k)=0.
         uz(i,1,k)=0.+dpdzy1(i,k)
      enddo
      enddo
      do k=1,xsize(3)
      do i=1,xsize(1)
         ux(i,xsize(2),k)=0.+dpdxyn(i,k)
         uy(i,xsize(2),k)=0.
         uz(i,xsize(2),k)=0.+dpdzyn(i,k)
      enddo
      enddo
   else
!find j=1 and j=ny
      if (xstart(2)==1) then
         do k=1,xsize(3)
         do i=1,xsize(1)
            ux(i,1,k)=0.+dpdxy1(i,k)
            uy(i,1,k)=0.
            uz(i,1,k)=0.+dpdzy1(i,k)
         enddo
         enddo
      endif
!      print *,nrank,xstart(2),ny-(nym/p_row)
       if (ny-(nym/dims(1))==xstart(2)) then
         do k=1,xsize(3)
         do i=1,xsize(1)
            ux(i,xsize(2),k)=0.+dpdxyn(i,k)
            uy(i,xsize(2),k)=0.
            uz(i,xsize(2),k)=0.+dpdzyn(i,k)
         enddo
         enddo
      endif
 
   endif
   endif
!********itype == 5******JICF*******************************   
   if (itype.eq.5) then
   ! determine the processor grid in use
   !call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, &
   !      dims, dummy_periods, dummy_coords, code)
  ujicf=0.
  njicf=0.
  temp=0.
!********bottom wall*************************************
  if (v_jicf==0.) then !nojet
    if (xstart(2)==1) then
      do k=1,xsize(3)
      do i=1,xsize(1)
         ux(i,1,k)=0.+dpdxy1(i,k)
         uy(i,1,k)=0.
         uz(i,1,k)=0.+dpdzy1(i,k)
      enddo
      enddo
    endif
  else  !jet
    njicf1=0.
    ujicf1=0.
    temp1=0.
    if (xstart(2)==1) then
      do k=1,xsize(3)
      do i=1,xsize(1)
	ux(i,1,k)=0.+dpdxy1(i,k)
        !check the hole
        tempa=xjicf-(i+xstart(1)-1)
        tempb=nz/2-(k+xstart(3)-1)
        r2=dx2*real(tempa*tempa)+dz2*real(tempb*tempb)
        if (r2 .le. rjicf) then
	  r2=2.*v_jicf*(1.-r2/rjicf)
	  !uy(i,1,k)=min(r2,(itime-ifirst)*r2/10000)
          uy(i,1,k)=r2
          temp1=max(temp1,uy(i,1,k))
          ujicf1=ujicf1+uy(i,1,k)*dx*dz
          njicf1=njicf1+1
        else
          uy(i,1,k)=0.
        endif
        uz(i,1,k)=0.+dpdzy1(i,k)
      enddo
      enddo
    endif
    call MPI_ALLREDUCE(ujicf1,ujicf,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(njicf1,njicf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,code)
    call MPI_REDUCE(temp1,temp,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
    if ((nrank==0) .and. (print_flag==1)) print *,"njicf",njicf,"ujicf",ujicf,"max V",temp
  endif  !  if (v_jicf==0) then
!      print *,nrank,xstart(2),ny-(nym/p_row)

!********top wall*************************************
  if (xend(2)==ny) then
    do k=1,xsize(3)
    do i=1,xsize(1)
      ux(i,xsize(2),k)=0.+dpdxyn(i,k)
      uy(i,xsize(2),k)=0.
      uz(i,xsize(2),k)=0.+dpdzyn(i,k)
    enddo
    enddo
  endif
endif  !itype == 5
endif  !NCLY==2

!      print *,nrank,'dims1',dims(1),'dims2',dims(2)
!      print *,nrank, njicf
!****************************************************
!********NCLZ==2*************************************
!****************************************************
!****************************************************

!##################################################### 

return
end subroutine pre_correc
