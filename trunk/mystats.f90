module user_stats

use decomp_2d, only : mytype, DECOMP_INFO

implicit none

interface update_user_stats
   module procedure update_user_stats_non_phi
   module procedure update_user_stats_oui_phi
end interface update_user_stats

logical, save :: bool_user_stat
integer, save :: beg_stat
TYPE(DECOMP_INFO) :: decomp_user_stats
integer, save :: xst1, xst2, xst3
integer, save :: xen1, xen2, xen3

! Ordre 1
real(mytype), save, allocatable, dimension(:,:,:) :: um,vm,wm

! Ordre 2
real(mytype), save, allocatable, dimension(:,:,:) :: uum,vvm,wwm
real(mytype), save, allocatable, dimension(:,:,:) :: uvm,uwm,vwm

! Scalaire
! Ordre 1
real(mytype), save, allocatable, dimension(:,:,:) :: phim
! Ordre 2
real(mytype), save, allocatable, dimension(:,:,:) :: uphim, vphim, wphim, phiphim

contains

subroutine allocate_user_stats(nx_user, ny_user, nz_user, phG,ph1,ph2,ph3,ph4)

use decomp_2d, only : decomp_info_init, DECOMP_INFO
use param, only : iscalar

implicit none

integer, intent(in) :: nx_user, ny_user, nz_user
TYPE(DECOMP_INFO), intent(in) :: phG,ph1,ph2,ph3,ph4

call decomp_info_init(nx_user, ny_user, nz_user, decomp_user_stats)

#include "user_module_diff.f90"

allocate(um     (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(vm     (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(wm     (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(uum    (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(vvm    (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(wwm    (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(uvm    (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(uwm    (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(vwm    (xst1:xen1,xst2:xen2,xst3:xen3))

um=0.; vm=0.; wm=0.;
uum=0.; vvm=0.; wwm=0.;
uvm=0.; uwm=0.; vwm=0.;

! Scalaire
if (iscalar.eq.1) then
  ! Ordre 1
  allocate(phim(xst1:xen1,xst2:xen2,xst3:xen3))
  ! Ordre 2
  allocate(uphim(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(vphim(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(wphim(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(phiphim(xst1:xen1,xst2:xen2,xst3:xen3))
  !
  phim=0.;
  uphim=0.; vphim=0.; wphim=0.; phiphim=0.;
  !
endif

return
end subroutine allocate_user_stats

subroutine update_user_stats_non_phi(ux,uy,uz,phG,ph1,ph2,ph3,ph4)

use decomp_2d, only : xsize, nrank, DECOMP_INFO
use param, only : itime

implicit none

TYPE(DECOMP_INFO), intent(in) :: phG,ph1,ph2,ph3,ph4
real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ux,uy,uz

! Local variable
integer :: coeff

coeff=itime-beg_stat
if (coeff.le.0) then
  if (nrank.eq.0) print*,'ERROR : cannot start to compute statistics at ',itime,'. Wait till ',beg_stat
  return
endif

um = ((coeff-1)*1./coeff)*um + (1./coeff)*fine_interpol(ux)
vm = ((coeff-1)*1./coeff)*vm + (1./coeff)*fine_interpol(uy)
wm = ((coeff-1)*1./coeff)*wm + (1./coeff)*fine_interpol(uz)

uum = ((coeff-1)*1./coeff)*uum + (1./coeff)*fine_interpol(ux*ux)
vvm = ((coeff-1)*1./coeff)*vvm + (1./coeff)*fine_interpol(uy*uy)
wwm = ((coeff-1)*1./coeff)*wwm + (1./coeff)*fine_interpol(uz*uz)

uvm = ((coeff-1)*1./coeff)*uvm + (1./coeff)*fine_interpol(ux*uy)
uwm = ((coeff-1)*1./coeff)*uwm + (1./coeff)*fine_interpol(ux*uz)
vwm = ((coeff-1)*1./coeff)*vwm + (1./coeff)*fine_interpol(uy*uz)

return
end subroutine update_user_stats_non_phi

subroutine update_user_stats_oui_phi(ux,uy,uz,phi,phG,ph1,ph2,ph3,ph4)

use decomp_2d, only : xsize, nrank, DECOMP_INFO
use param, only : itime

implicit none

TYPE(DECOMP_INFO), intent(in) :: phG,ph1,ph2,ph3,ph4
real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ux,uy,uz
real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: phi

! Local variable
integer :: coeff

coeff=itime-beg_stat
if (coeff.le.0) then
  if (nrank.eq.0) print*,'ERROR : cannot start to compute statistics at ',itime,'. Wait till ',beg_stat
  return
endif

call update_user_stats(ux,uy,uz,phG,ph1,ph2,ph3,ph4)

phim = ((coeff-1)*1./coeff)*phim + (1./coeff)*fine_interpol(phi)

uphim = ((coeff-1)*1./coeff)*uphim + (1./coeff)*fine_interpol(ux*phi)
vphim = ((coeff-1)*1./coeff)*vphim + (1./coeff)*fine_interpol(uy*phi)
wphim = ((coeff-1)*1./coeff)*wphim + (1./coeff)*fine_interpol(uz*phi)
phiphim = ((coeff-1)*1./coeff)*phiphim + (1./coeff)*fine_interpol(phi*phi)

return
end subroutine update_user_stats_oui_phi

subroutine pre_update_user_stats(phG,ph1,ph2,ph3,ph4)

  use decomp_2d, only : xsize, ysize, zsize, DECOMP_INFO
  use param, only : istret, iscalar, iimplicit
  use variables
  use var

  implicit none

  TYPE(DECOMP_INFO), intent(in) :: phG,ph1,ph2,ph3,ph4
  logical :: myflag
  integer :: i,j,k
  real(mytype), allocatable, dimension(:,:,:) :: tmp_dphidx,tmp_dphidy,tmp_dphidz,tmp_deltaphi
  ! Pressure interpolation
  real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: tmp_pres
  real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),zsize(3)) :: myta3,mydi3
  real(mytype),dimension(ph3%yst(1):ph3%yen(1),nymsize,ysize(3)) :: myta2
  real(mytype),dimension(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)) :: mytb2,mydi2
  real(mytype),dimension(nxmsize,xsize(2),xsize(3)) :: myta1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: mytb1

  
if (.false.) then
  ! Compute laplacian temperature (budget uPhi, vPhi, wPhi) and grad(phi)
  if (iscalar.eq.1) then

    ! Compute laplacian phi
    allocate(tmp_deltaphi(xsize(1),xsize(2),xsize(3)))
    ! Compute grad phi
    allocate(tmp_dphidx(xsize(1),xsize(2),xsize(3)))
    allocate(tmp_dphidy(xsize(1),xsize(2),xsize(3)))
    allocate(tmp_dphidz(xsize(1),xsize(2),xsize(3)))
    !
    call derxxt(ta1,phi1,di1,sx,sfxt,ssxt,swxt,xsize(1),xsize(2),xsize(3),0)
    call derx(tmp_dphidx,phi1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
    !
    call transpose_x_to_y(phi1,phi2)
    !
    myflag=.false.
    if (iimplicit.eq.1) myflag=.true.
    if (myflag) iimplicit=0
    call deryyt(ta2,phi2,di2,sy,sfyt,ssyt,swyt,ysize(1),ysize(2),ysize(3),0)
    if (myflag) iimplicit=1
    call dery (tc2,phi2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
    if (istret.ne.0) then 
      do k=1,ysize(3)
      do j=1,ysize(2)
      do i=1,ysize(1)
        ta2(i,j,k)=ta2(i,j,k)*pp2y(j)-pp4y(j)*tc2(i,j,k) ! dyy dans ta2
      enddo
      enddo
      enddo
    endif
    call transpose_y_to_x(tc2,tmp_dphidy) ! dTdy
    !
    call transpose_y_to_z(phi2,phi3)
    !
    call derzzt (ta3,phi3,di3,sz,sfzt,sszt,swzt,zsize(1),zsize(2),zsize(3),0)
    call transpose_z_to_y(ta3,tc2) ! dzz dans tc2
    call derz (ta3,phi3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
    call transpose_z_to_y(ta3,tb2) ! dTdz dans tb2
    call transpose_y_to_x(tb2,tmp_dphidz) ! dTdz
    !
    tc2=tc2+ta2 ! dzz + dyy
    call transpose_y_to_x(tc2,tc1) ! dzz + dyy
    tmp_deltaphi=ta1+tc1 ! dzz + dyy + dxx
    !
  endif

  ! Compute velocity gradient
  call derx (ta1,ux1,di1,sx,ffx ,fsx ,fwx ,xsize(1),xsize(2),xsize(3),0)
  call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
  call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
  ! ta1=dudx / tb1=dvdx / tc1=dwdx
  call transpose_x_to_y(ux1,ux2)
  call transpose_x_to_y(uy1,uy2)
  call transpose_x_to_y(uz1,uz2)
  call dery (ta2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
  call dery (tb2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
  call dery (tc2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
! ta2=dudy / tb2=dvdy / tc2=dwdy
  call transpose_y_to_z(ux2,ux3)
  call transpose_y_to_z(uy2,uy3)
  call transpose_y_to_z(uz2,uz3)
  call derz (ta3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
  call derz (tb3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
  call derz (tc3,uz3,di3,sz,ffz ,fsz ,fwz ,zsize(1),zsize(2),zsize(3),0)
  ! ta3=dudz / tb3=dvdz / tc3=dwdz
  ! Back in y
  call transpose_z_to_y(ta3,td2)!dudz
  call transpose_z_to_y(tb3,te2)!dvdz
  call transpose_z_to_y(tc3,tf2)!dwdz
  ! Back to x
  call transpose_y_to_x(ta2,td1)!dudy
  call transpose_y_to_x(tb2,te1)!dvdy
  call transpose_y_to_x(tc2,tf1)!dwdy
  call transpose_y_to_x(td2,tg1)!dudz
  call transpose_y_to_x(te2,th1)!dvdz
  call transpose_y_to_x(tf2,ti1)!dwdz

  if (.true.) then
    !WORK Z-PENCILS
    call interiz6(myta3,pp3,mydi3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
         (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
    !WORK Y-PENCILS
    call transpose_z_to_y(myta3,myta2,ph3) !nxm nym nz
    call interiy6(mytb2,myta2,mydi2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
         (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
    !WORK X-PENCILS
    call transpose_y_to_x(mytb2,myta1,ph2) !nxm ny nz
    call interi6(tmp_pres,myta1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
         nxmsize,xsize(1),xsize(2),xsize(3),1)
    !The pressure field on the main mesh is in tmp_pres
  endif
endif
  
  if (iscalar.eq.1) then
    call update_user_stats(ux1,uy1,uz1,phi1,phG,ph1,ph2,ph3,ph4)
  else
    call update_user_stats(ux1,uy1,uz1,phG,ph1,ph2,ph3,ph4)
  endif

return
end subroutine pre_update_user_stats

subroutine write_user_stats(phG,ph1,ph2,ph3,ph4)

use decomp_2d, only : DECOMP_INFO
use decomp_2d_io, only : decomp_2d_write_plane, decomp_2d_write_one
use param, only : iscalar

implicit none

TYPE(DECOMP_INFO), intent(in) :: phG,ph1,ph2,ph3,ph4

call decomp_2d_write_one(1,um,'um.dat',decomp_user_stats)
call decomp_2d_write_one(1,vm,'vm.dat',decomp_user_stats)
call decomp_2d_write_one(1,wm,'wm.dat',decomp_user_stats)

call decomp_2d_write_one(1,uum,'uum.dat',decomp_user_stats)
call decomp_2d_write_one(1,vvm,'vvm.dat',decomp_user_stats)
call decomp_2d_write_one(1,wwm,'wwm.dat',decomp_user_stats)

call decomp_2d_write_one(1,uvm,'uvm.dat',decomp_user_stats)
call decomp_2d_write_one(1,uwm,'uwm.dat',decomp_user_stats)
call decomp_2d_write_one(1,vwm,'vwm.dat',decomp_user_stats)

if (iscalar.eq.1) then
  call decomp_2d_write_one(1,phim,'phim.dat',decomp_user_stats)
  !
  call decomp_2d_write_one(1,uphim,'uphim.dat',decomp_user_stats)
  call decomp_2d_write_one(1,vphim,'vphim.dat',decomp_user_stats)
  call decomp_2d_write_one(1,wphim,'wphim.dat',decomp_user_stats)
  call decomp_2d_write_one(1,phiphim,'phiphim.dat',decomp_user_stats)
  !
endif

return
end subroutine write_user_stats

function fine_interpol(var)

  USE decomp_2d, only : mytype, xsize

  implicit none

  real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: var
  real(mytype), dimension(xst1:xen1,xst2:xen2,xst3:xen3) :: fine_interpol
  integer :: i,j,k

  fine_interpol=0.

  do k=xst3,xen3
  do j=xst2,xen2
  do i=xst1,xen1
    fine_interpol(i,j,k)=var( (i-xst1+1)*xsize(1)/(xen1-xst1+1), &
			      (j-xst2+1)*xsize(2)/(xen2-xst2+1), &
			      (k-xst3+1)*xsize(3)/(xen3-xst3+1) )
  enddo
  enddo
  enddo

end function fine_interpol

subroutine read_user_stats(phG,ph1,ph2,ph3,ph4)
  !
  ! Read average when using a restart
  !
  use decomp_2d, only : DECOMP_INFO
  use decomp_2d_io, only : decomp_2d_read_one
  use param, only : iscalar

  implicit none

  TYPE(DECOMP_INFO), intent(in) :: phG,ph1,ph2,ph3,ph4

  character(len=30) :: filename

  filename='um.dat'
  call read_one_user_stat(filename,um)
  filename='vm.dat'
  call read_one_user_stat(filename,vm)
  filename='wm.dat'
  call read_one_user_stat(filename,wm)

  filename='uum.dat'
  call read_one_user_stat(filename,uum)
  filename='vvm.dat'
  call read_one_user_stat(filename,vvm)
  filename='wwm.dat'
  call read_one_user_stat(filename,wwm)

  filename='uvm.dat'
  call read_one_user_stat(filename,uvm)
  filename='uwm.dat'
  call read_one_user_stat(filename,uwm)
  filename='vwm.dat'
  call read_one_user_stat(filename,vwm)

  if (iscalar.eq.1) then
    filename='phim.dat'
    call read_one_user_stat(filename,phim)
    !
    filename='uphim.dat'
    call read_one_user_stat(filename,uphim)
    filename='vphim.dat'
    call read_one_user_stat(filename,vphim)
    filename='wphim.dat'
    call read_one_user_stat(filename,wphim)
    filename='phiphim.dat'
    call read_one_user_stat(filename,phiphim)
    !
  endif

return
end subroutine read_user_stats

subroutine read_one_user_stat(filename,var)

  USE decomp_2d, only : mytype, xsize, nrank
  use decomp_2d_io, only : decomp_2d_read_one      !for 3-D statistic
  !use decomp_2d_io, only : decomp_2d_read_plane    !for plane statistic

  use mpi

  implicit none

  real(mytype), dimension(xst1:xen1,xst2:xen2,xst3:xen3) :: var
  character(len=30) :: filename

  logical file_exist
  real(mytype), dimension(xst1:xen1,xst2:xen2,xst3:xen3) :: tmp	!for 3-D statistic
  !real(mytype), dimension(xst1:xst1,xst2:xen2,xst3:xen3) :: tmp		!for plane statistic
  integer i,j,k

  INQUIRE(FILE=filename, EXIST=file_exist)

  if (file_exist) then
    !
    call decomp_2d_read_one(1,tmp,filename,decomp_user_stats)		!for 3-D statistic
    !call decomp_2d_read_plane(1,tmp,filename,decomp_user_stats)		!for plane statistic
    !
    do k=xst3,xen3
      do j=xst2,xen2
        do i=xst1,xen1
          var(i,j,k)=tmp(i,j,k)		!for 3-D statistic
          !var(i,j,k)=tmp(xst1,j,k)		!for plane statistic
        end do
      end do
    end do
    !
    if (nrank.eq.0) print*,'user_stats : reading average for file ',trim(filename),' done'
    !
  else
    !
    var = 0.
    !
    if (nrank.eq.0) print*,'user_stats : average file ',trim(filename),' is not existing, initial value is set to 0'
    !
  endif

  return
end subroutine read_one_user_stat

subroutine plot_my_stats()

use decomp_2d
use param
use variables
use var
USE MPI

implicit none

integer :: j,i,k,code
real(mytype) :: u_tau1,Re_tau1,sum1,sum

call transpose_x_to_y(ux1,ux2)
call dery (ta2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)

u_tau1=0
sum1=0
Re_tau1=0
sum=0

if (itype==5 .and. rezone>0) then
  do k=1,ysize(3)
  do i=1,ysize(1)
    if(ystart(1)+i-1 .le. rezone) then
      u_tau1=u_tau1+abs(ta2(i,1,k))+abs(ta2(i,ny,k))
      sum1=sum1+ta2(i,1,k)+ta2(i,ny,k)
    endif
  enddo
  enddo

  call MPI_ALLREDUCE(u_tau1,Re_tau1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
  call MPI_ALLREDUCE(sum1,sum,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
  Re_tau1=sqrt(Re_tau1/rezone/nz/2/xnu)
  sum=sum/rezone/nz/2
  if (print_flag==1 .and. nrank==0) print *,'Re_tau',Re_tau1,'top-bottom diff',sum
  
  Re_tau=(Re_tau*(itime-ifirst)+Re_tau1)/(itime-ifirst+1)
  if (print_flag==1 .and. nrank==0) print *,'Time-averaged h+',Re_tau
else if (itype==2) then
  do k=1,ysize(3)
  do i=1,ysize(1)
    u_tau1=u_tau1+abs(ta2(i,1,k))+abs(ta2(i,ny,k))
    sum1=sum1+ta2(i,1,k)+ta2(i,ny,k)
  enddo
  enddo

  call MPI_ALLREDUCE(u_tau1,Re_tau1,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
  call MPI_ALLREDUCE(sum1,sum,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
  Re_tau1=sqrt(Re_tau1/rezone/nz/2/xnu)
  sum=sum/rezone/nz/2
  if (print_flag==1 .and. nrank==0) print *,'Re_tau',Re_tau1,'top-bottom diff',sum
  
  Re_tau=(Re_tau*(itime-ifirst)+Re_tau1)/(itime-ifirst+1)
  if (print_flag==1 .and. nrank==0) print *,'Time-averaged h+',Re_tau
endif

call transpose_y_to_x(ta2,td1)!dudy

return
end subroutine

end module user_stats