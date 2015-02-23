!#define my_mod_solide

!************************************************************
!
subroutine scalar_exp(ux1,uy1,uz1,phi1,phis1,phiss1,di1,ta1,tb1,tc1,td1,&
     uy2,uz2,phi2,di2,ta2,tb2,tc2,td2,uz3,phi3,di3,ta3,tb3,epsi)
!
!************************************************************

USE param
USE variables
USE decomp_2d

implicit none

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,phi1,phis1,&
                                              phiss1,di1,ta1,tb1,tc1,td1,epsi
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2,uz2,phi2,di2,ta2,tb2,tc2,td2,temp1,temp2,temp3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: uz3,phi3,di3,ta3,tb3
real(mytype),dimension(ysize(1),2*ysize(2),ysize(3)) :: tt2

integer :: ijk,nvect1,nvect2,nvect3,i,j,k,l,nxyz,tempa,tempb
real(mytype) :: x,y,z,r2
logical,dimension(ysize(1),ysize(3)) :: injet

nvect1=xsize(1)*xsize(2)*xsize(3)
nvect2=ysize(1)*ysize(2)*ysize(3)
nvect3=zsize(1)*zsize(2)*zsize(3)

!call transpose_x_to_y(ux1,temp1)
!call transpose_x_to_y(uy1,temp2)
!call transpose_x_to_y(uz1,temp3)
!if (print_flag==1 .and. ystart(1)==1 .and. ystart(3)==1) then
!  print *,'ux',temp1(1,:,1)
!  print *,'uy',temp2(1,:,1)
!  print *,'uz',temp3(1,:,1)
!endif

!X PENCILS
100 do ijk=1,nvect1
   ta1(ijk,1,1)=ux1(ijk,1,1)*phi1(ijk,1,1)
enddo
call derx (tb1,ta1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
call derxx (ta1,phi1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)

!call transpose_x_to_y(tb1,temp1)
!call transpose_x_to_y(ta1,temp2)
!if (print_flag==1 .and. ystart(1)==1 .and. ystart(3)==1) then
!  print *,'dTUdx',temp1(1,:,1)
!  print *,'dTdxx',temp2(1,:,1)
!endif

call transpose_x_to_y(phi1,phi2)
call transpose_x_to_y(uy1,uy2)
call transpose_x_to_y(uz1,uz2)


!Y PENCILS
do ijk=1,nvect2
  ta2(ijk,1,1)=uy2(ijk,1,1)*phi2(ijk,1,1)
enddo
call dery (tb2,ta2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)

if (istret.ne.0) then 
  call deryy (ta2,phi2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
  call dery (tc2,phi2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
  do k=1,ysize(3)
  do j=1,ysize(2)
  do i=1,ysize(1)
    ta2(i,j,k)=ta2(i,j,k)*pp2y(j)-pp4y(j)*tc2(i,j,k)
  enddo
  enddo
  enddo
else
  call deryy (ta2,phi2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1) 
endif

!if (print_flag==1 .and. ystart(1)==1 .and. ystart(3)==1) then
!  print *,'dTUdy',tb2(1,:,1)
!  print *,'dTdyy',ta2(1,:,1)
!endif

call transpose_y_to_z(phi2,phi3)
call transpose_y_to_z(uz2,uz3)

!Z PENCILS
do ijk=1,nvect3
   ta3(ijk,1,1)=uz3(ijk,1,1)*phi3(ijk,1,1)
enddo
call derz (tb3,ta3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
call derzz (ta3,phi3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)

call transpose_z_to_y(ta3,tc2)
call transpose_z_to_y(tb3,td2)

!if (print_flag==1 .and. ystart(1)==1 .and. ystart(3)==1) then
!   print *,'dTUdz',td2(1,:,1)
!   print *,'dTdzz',tc2(1,:,1)
!endif

!Y PENCILS ADD TERMS
do ijk=1,nvect2
   tc2(ijk,1,1)=tc2(ijk,1,1)+ta2(ijk,1,1)
   td2(ijk,1,1)=td2(ijk,1,1)+tb2(ijk,1,1)
enddo

call transpose_y_to_x(tc2,tc1)
call transpose_y_to_x(td2,td1)

!X PENCILS ADD TERMS
do ijk=1,nvect1
   ta1(ijk,1,1)=ta1(ijk,1,1)+tc1(ijk,1,1) !SECOND DERIVATIVE
   tb1(ijk,1,1)=tb1(ijk,1,1)+td1(ijk,1,1) !FIRST DERIVATIVE
enddo

! call transpose_x_to_y(tb1,temp1)
! call transpose_x_to_y(ta1,temp2)
! if (print_flag==1 .and. ystart(1)==1 .and. ystart(3)==1) then
!   print *,'FIRST DERIVATIVE',temp1(1,:,1)
!   print *,'SECOND DERIVATIVE',temp2(1,:,1)
! endif

do ijk=1,nvect1
   ta1(ijk,1,1)=xnu/sc*ta1(ijk,1,1)-tb1(ijk,1,1)-xnu/sc*ux1(ijk,1,1)
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Kasagi Source term !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!
!ta1=ta1-xnu/sc*ux1

do ijk=1,nvect1
   ta1(ijk,1,1)=(1.-epsi(ijk,1,1))*ta1(ijk,1,1)
enddo

!TIME ADVANCEMENT
nxyz=xsize(1)*xsize(2)*xsize(3)  

if ((nscheme.eq.1).or.(nscheme.eq.2)) then
  if ((nscheme.eq.1.and.itime.eq.1.and.ilit.eq.0).or.&
	(nscheme.eq.2.and.itr.eq.1)) then
    do ijk=1,nxyz
      if (epsi(ijk,1,1)==0.) then
	phi1(ijk,1,1)=gdt(itr)*ta1(ijk,1,1)+phi1(ijk,1,1)
      else 
	phi1(ijk,1,1)=0.
      endif
      phis1(ijk,1,1)=ta1(ijk,1,1)          
    enddo
  else
    do ijk=1,nxyz
      if (epsi(ijk,1,1)==0.) then      
	phi1(ijk,1,1)=adt(itr)*ta1(ijk,1,1)+bdt(itr)*phis1(ijk,1,1)+phi1(ijk,1,1)
      else
	phi1(ijk,1,1)=0.
      endif
      phis1(ijk,1,1)=ta1(ijk,1,1)          
    enddo
  endif
endif

if (nscheme.eq.3) then 
if (nrank==0) print *,'Not ready'
stop 
endif

if (nscheme==4) then
   if ((itime.eq.1).and.(ilit.eq.0)) then
      if (nrank==0) print *,'start with Euler',itime
      do ijk=1,nxyz !start with Euler
         phi1(ijk,1,1)=dt*ta1(ijk,1,1)+phi1(ijk,1,1)
         phis1(ijk,1,1)=ta1(ijk,1,1)          
      enddo
   else
      if  ((itime.eq.2).and.(ilit.eq.0)) then
         if (nrank==0) print *,'then with AB2',itime
         do ijk=1,nxyz
            phi1(ijk,1,1)=1.5*dt*ta1(ijk,1,1)-0.5*dt*phis1(ijk,1,1)+phi1(ijk,1,1)
            phiss1(ijk,1,1)=phis1(ijk,1,1)
            phis1(ijk,1,1)=ta1(ijk,1,1)
         enddo 
      else
         do ijk=1,nxyz
            phi1(ijk,1,1)=adt(itr)*ta1(ijk,1,1)+bdt(itr)*phis1(ijk,1,1)+&
                 cdt(itr)*phiss1(ijk,1,1)+phi1(ijk,1,1)
            phiss1(ijk,1,1)=phis1(ijk,1,1)
            phis1(ijk,1,1)=ta1(ijk,1,1)
         enddo
      endif
   endif
endif

call transpose_x_to_y(phi1,phi2)
tt2=0.
injet=.false.

if (ncly==2) then
if (itype.eq.2) then
!********bottom wall*************************************
  phi2(:,1,:)=g_0
!********top wall*************************************
  phi2(:,ny,:)=g_n
else if (itype.eq.5) then
!********bottom wall*************************************
  if (v_jicf==0.) then
    phi2(:,1,:)=0. !g_0
  else
    do k=1,ysize(3)
    do i=1,ysize(1)
      !check the hole
      tempa=xjicf-(i+ystart(1)-1)
      tempb=nz/2-(k+ystart(3)-1)
      r2=dx2*real(tempa*tempa)+dz2*real(tempb*tempb)
      if (r2 .le. rjicf) then
        phi2(i,1,k)=1.
        injet(i,k)=.true.
      else
        phi2(i,1,k)=0. !g_0
      endif
    enddo
    enddo
  endif
!********top wall*************************************
  phi2(:,ny,:)=0. !g_n
endif
endif

do k=1,ysize(3)
do l=1,ny
do j=1,ny+ny
do i=1,ysize(1)
  if (injet(i,k)) then
    tt2(i,j,k)=tt2(i,j,k)+scalar_left_jet(j,l)*phi2(i,l,k)
  else
    tt2(i,j,k)=tt2(i,j,k)+scalar_left_main(j,l)*phi2(i,l,k)
  endif
enddo
enddo
enddo
enddo

do k=1,ysize(3)
do j=1,ysize(2)
do i=1,ysize(1)
    phi2(i,j,k)=tt2(i,j,k)
enddo
enddo
enddo
endif
! if (print_flag==1 .and. ystart(1)==1 .and. ystart(3)==1) then
!   !print *,'dTdy',tt2(1,ny+1:ny*2:4,1)
!   print *,'T',phi2(1,1:ny:8,1)
! endif


call transpose_y_to_x(phi2,phi1)

return
end subroutine scalar_exp
 
 
!********************************************************************
!
subroutine scalar_schemes_exp(alpha_00,beta_00,alpha_nn,beta_nn,scalar_left) 
!
!********************************************************************
USE derivY
USE variables
USE param
USE MPI,only : MPI_WTIME
USE decomp_2d, only : nrank

implicit none

integer :: i,j,ipiv(2*ny),info

real(mytype) :: t1,t2
real(mytype), intent(in) :: alpha_00,beta_00,alpha_nn,beta_nn
real(mytype), dimension(2*ny,2*ny) :: temp
real(mytype), dimension(2*ny,2*ny), intent(out) :: scalar_left

  t1 = MPI_WTIME()
  scalar_left=0.
  temp=0.
  do i=1,2*ny
     temp(i,i)=1.
  enddo
  scalar_left(1,1   )= alpha_00
  scalar_left(1,ny+1)=-beta_00
  do i=2,ny-1
    scalar_left(i,i)=1.
  enddo
  scalar_left(ny,ny   )=alpha_nn
  scalar_left(ny,ny+ny)=beta_nn
  
  scalar_left(ny+1,1)=-af1y
  scalar_left(ny+1,2)=-bf1y
  scalar_left(ny+1,3)=-cf1y
  scalar_left(ny+2,1)= af2y
  scalar_left(ny+2,3)=-af2y
  do j=3,ny-2
    scalar_left(ny+j,j-2)= bfjy
    scalar_left(ny+j,j-1)= afjy
    scalar_left(ny+j,j+1)=-afjy
    scalar_left(ny+j,j+2)=-bfjy
  enddo
  scalar_left(ny+ny-1,ny-2)= afmy
  scalar_left(ny+ny-1,ny  )=-afmy
  scalar_left(ny+ny  ,ny-2)= cfny
  scalar_left(ny+ny  ,ny-1)= bfny
  scalar_left(ny+ny  ,ny  )= afny

  if (istret.ne.0) then
    scalar_left(ny+1,ny+1)=fcy(1)/ppy(1)
    scalar_left(ny+1,ny+2)=ffy(1)/ppy(2)
    do j=2,ny-1
      scalar_left(ny+j,ny+j-1)=fby(j-1)/ppy(j-1)
      scalar_left(ny+j,ny+j  )=fcy(j)/ppy(j)
      scalar_left(ny+j,ny+j+1)=ffy(j)/ppy(j+1)
    enddo
    scalar_left(ny+ny,ny+ny-1)=fby(ny-1)/ppy(ny-1)
    scalar_left(ny+ny,ny+ny)=fcy(ny)/ppy(ny)
  else
    scalar_left(ny+1,ny+1)=fcy(1)
    scalar_left(ny+1,ny+2)=ffy(1)
    do j=2,ny-1
      scalar_left(ny+j,ny+j-1)=fby(j-1)
      scalar_left(ny+j,ny+j  )=fcy(j)
      scalar_left(ny+j,ny+j+1)=ffy(j)
    enddo
    scalar_left(ny+ny,ny+ny-1)=fby(ny-1)
    scalar_left(ny+ny,ny+ny)=fcy(ny)
  endif

  t2=MPI_WTIME()-t1
  if (nrank==0) print *,'   - Matrix constructed ',t2,' seconds'
!print *,scalar_left
  t1 = MPI_WTIME()

#ifdef DOUBLE_PREC
  call dgesv(2*ny,2*ny,scalar_left,2*ny,ipiv,temp,2*ny,info)
#else
  call sgesv(2*ny,2*ny,scalar_left,2*ny,ipiv,temp,2*ny,info)
#endif

  scalar_left=temp
  t2=MPI_WTIME()-t1
  if (nrank==0) print *,'   - Matrix inversed ',t2,' seconds'

  return
end subroutine scalar_schemes_exp
