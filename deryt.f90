!********************************************************************
!
subroutine dery_T(ty,uy,ry,sy,ffy,fsy,fwy,ppy,nx,ny,nz,npaire) 
!
!********************************************************************
!!For Neumann or Robin Boundary conditions
!!
!! Robin boundary condition on temperature
!! alpha * T + beta * dT/dn = g
!! alpha=1, beta=0 is dirichlet
!! alpha=0, beta=1 is neumann
!!
!! WARNING ATTENTION ACHTUNG WARNING ATTENTION ACHTUNG
!!
!! beta is the coefficient for NORMAL derivative :
!!
!! alpha_0*T(0) - beta_0*dTdy(0)=g_0
!! alpha_n*T(L) + beta_n*dTdy(L)=g_n
!!
!! WARNING ATTENTION ACHTUNG WARNING ATTENTION ACHTUNG
!!


USE param
USE derivY
use decomp_2d, only : nrank

implicit none

integer :: nx,ny,nz,i,j,k,npaire
real(mytype), dimension(nx,ny,nz) :: ty,uy
real(mytype), dimension(nx,ny,nz) :: ry
real(mytype), dimension(nx,nz)  :: sy
real(mytype), dimension(ny) :: ffy,fsy,fwy,ppy

if (ncly.ne.2) then
  if (nrank==0) print *,'dery_T for ncly==0/1 Not ready '
  stop
endif

ty(:,1,:)=af1y*uy(:,1,:)+bf1y*uy(:,2,:)+cf1y*uy(:,3,:)
ty(:,2,:)=af2y*(uy(:,3,:)-uy(:,1,:))
do j=3,ny-2
  ty(:,j,:)=afjy*(uy(:,j+1,:)-uy(:,j-1,:))&
    +bfjy*(uy(:,j+2,:)-uy(:,j-2,:))
enddo
ty(:,ny-1,:)=afmy*(uy(:,ny,:)-uy(:,ny-2,:))
ty(:,ny,:)=-afny*uy(:,ny,:)-bfny*uy(:,ny-1,:)-cfny*uy(:,ny-2,:)

if (beta_0 .ne. 0.) ty(:,1,:)=g_0-alpha_0*uy(:,1,:)
if (beta_n .ne. 0.) ty(:,ny,:)=g_n-alpha_n*uy(:,ny,:)

do j=2,ny
  ty(:,j,:)=ty(:,j,:)-ty(:,j-1,:)*fsy(j)
enddo
ty(:,ny,:)=ty(:,ny,:)*fwy(ny)
do j=ny-1,1,-1
  ty(:,j,:)=(ty(:,j,:)-ffy(j)*ty(:,j+1,:))*fwy(j)
enddo

if (istret.ne.0) then
   do k=1,nz
   do j=1,ny
   do i=1,nx
      ty(i,j,k)=ty(i,j,k)*ppy(j)
   enddo
   enddo
   enddo
endif

return
end subroutine dery_T