module conjugate_ht
!
! Module for conjugate heat transfer
!
use decomp_2d, only : mytype, DECOMP_INFO

type(decomp_info), save :: mydecomp_bot, mydecomp_top
logical, save :: bool_conjugate_ht, bool_sol_stats
integer, save :: beg_stat_sol
integer, save :: ny_sol_bot, ny_sol_top
real(mytype), parameter :: sol_imp_var=0. ! 1=explicit, 0=implicit, 0.5=Crank-Nicolson
real(mytype), save :: ly_sol_bot, ly_sol_top
real(mytype), save :: repr_sol_bot, repr_sol_top
real(mytype), save :: fluxratio_bot, fluxratio_top
real(mytype), save, allocatable, dimension(:) :: yp_bot, yp_top
real(mytype), save, allocatable, dimension(:,:) :: val_bot, val_top
real(mytype), save, allocatable, dimension(:,:) :: flux_bot, flux_top
real(mytype), save, allocatable, dimension(:,:) :: d2pn_bot, d2pn_top
real(mytype), save, allocatable, dimension(:,:) :: invdiffy_bot, invdiffy_top
real(mytype), save, allocatable, dimension(:,:,:) :: cl_bot, cl_top
real(mytype), save, allocatable, dimension(:,:,:) :: temp_bot, temp_top
real(mytype), save, allocatable, dimension(:,:,:) :: bot_dm1, top_dm1
real(mytype), save, allocatable, dimension(:,:,:) :: bot_dm2, top_dm2
real(mytype), save, allocatable, dimension(:,:,:) :: bot_dm3, top_dm3
!
! stats
!
real(mytype), save, allocatable, dimension(:,:,:) :: tempm_1_bot, tempm_2_bot
real(mytype), save, allocatable, dimension(:,:,:) :: tempm_1_top, tempm_2_top

contains

subroutine conjugate_ht_init()
!
! Subroutine to initialize the module
!
  use decomp_2d, only : ysize, nrank, ystart
  use param, only : yly, ilit, twopi, dz, dx, dt

  implicit none

  ! Local variables
  integer :: i,j,k
  real(mytype) :: ytmp
  real(mytype) :: x, z, mysx, mysy, mysz
  real(mytype), dimension(ysize(1),ny_sol_bot+1+2,ysize(3)) :: ydiff_bot
  real(mytype), dimension(ysize(1),ny_sol_top+1+2,ysize(3)) :: ydiff_top

  if (nrank.eq.0) print *,'  User module for conjugate heat transfer :'
  if (nrank.eq.0) print *,'    - ',ny_sol_bot,' & ',ny_sol_top,' Chebyshev nodes top & bottom'

  ! ny_bot et ny_top : ordre de l'interpolation
  ! ny_bot+1 et ny_top+1 : nombres de points de chebyshev
  ! +2 points pour les bords
  !
  allocate(yp_bot(ny_sol_bot+1+2), yp_top(ny_sol_top+1+2))
  call cheb_nodes(yp_bot,ny_sol_bot,-ly_sol_bot,0.)
  call cheb_nodes(yp_top,ny_sol_top,yly,yly+ly_sol_top)
  if (nrank.eq.0) print *,'    - Chebyshev nodes computed'
  if (nrank.eq.0) then
   open(21,file='yp_bot.dat',form='formatted')
   do j=1,ny_sol_bot+1+2
    write(21,*) yp_bot(j)
   enddo
   close(21)
   open(21,file='yp_top.dat',form='formatted')
   do j=1,ny_sol_top+1+2
    write(21,*) yp_top(j)
   enddo
   close(21)
  endif
  !
  ! Compute coefficients for value at boundary
  !
  allocate(val_bot(2,ny_sol_bot+1),val_top(2,ny_sol_top+1))
  call cheb_edge(val_bot,yp_bot,ny_sol_bot)
  call cheb_edge(val_top,yp_top,ny_sol_top)
  if (nrank.eq.0) print *,'    - Boundary-value matrix computed'
  !
  ! Compute coefficients for flux at boundary
  !
  allocate(flux_bot(2,ny_sol_bot+1),flux_top(2,ny_sol_top+1))
  call cheb_flux(flux_bot,yp_bot,ny_sol_bot)
  call cheb_flux(flux_top,yp_top,ny_sol_top)
  if (nrank.eq.0) print *,'    - Boundary-flux matrix computed'
  !
  ! Compute coefficients of second derivative
  ! on Chebyshev nodes (y-direction)
  !
  allocate(d2pn_bot(ny_sol_bot+1,ny_sol_bot+1), d2pn_top(ny_sol_top+1,ny_sol_top+1))
  call cheb_der2(d2pn_bot,yp_bot,ny_sol_bot)
  call cheb_der2(d2pn_top,yp_top,ny_sol_top)
  if (nrank.eq.0) print *,'    - Y-diffusion matrix computed'
  if (nrank.eq.0) print *,'    - Computing inverse of matrix'
  allocate(invdiffy_bot(ny_sol_bot-1,ny_sol_bot-1))
  i = my_ydiff_mat_inv(d2pn_bot,ny_sol_bot,invdiffy_bot,flux_bot(1,:),flux_bot(2,:),repr_sol_bot)
  if (nrank.eq.0) print *,'      o Bottom matrix inverted : ', i
  allocate(invdiffy_top(ny_sol_top-1,ny_sol_top-1))
  i = my_ydiff_mat_inv(d2pn_top,ny_sol_top,invdiffy_top,flux_top(1,:),flux_top(2,:),repr_sol_top)
  if (nrank.eq.0) print *,'      o Top matrix inverted : ', i
  !
  allocate(cl_bot(ysize(1),1,ysize(3)))
  allocate(cl_top(ysize(1),1,ysize(3)))
  cl_bot=0.
  cl_top=0.
  !
  allocate(temp_bot(ysize(1),ny_sol_bot+1+2,ysize(3)))
  allocate(temp_top(ysize(1),ny_sol_top+1+2,ysize(3)))
  temp_bot=0.
  temp_top=0.
  allocate(bot_dm1(ysize(1),ny_sol_bot+1+2,ysize(3)))
  allocate(top_dm1(ysize(1),ny_sol_top+1+2,ysize(3)))
  bot_dm1=0.
  top_dm1=0.
  allocate(bot_dm2(ysize(1),ny_sol_bot+1+2,ysize(3)))
  allocate(top_dm2(ysize(1),ny_sol_top+1+2,ysize(3)))
  bot_dm2=0.
  top_dm2=0.
  allocate(bot_dm3(ysize(1),ny_sol_bot+1+2,ysize(3)))
  allocate(top_dm3(ysize(1),ny_sol_top+1+2,ysize(3)))
  bot_dm3=0.
  top_dm3=0.
  !
  ! 2DECOMP for bottom and top solid domain
  !
  call my2decomp_solide()
  !
  ! Initial value for temperature
  !
  if (ilit.eq.0) then
    !
    do k=1,ysize(3)
      z=(k-1+ystart(3)-1)*dz
      mysz=sin(twopi*z)
    do j=1,ny_sol_bot+1+2
      ytmp=yp_bot(j)
      mysy=sin(twopi*ytmp)
    do i=1,ysize(1)
      x=(i-1+ystart(1)-1)*dx
      mysx=sin(twopi*x)
      temp_bot(i,j,k)=ytmp*fluxratio_bot
    enddo
    enddo
    enddo
    !
    do k=1,ysize(3)
      z=(k-1+ystart(3)-1)*dz
      mysz=sin(twopi*z)
    do j=1,ny_sol_top+1+2
      ytmp=yp_top(j)
      mysy=sin(twopi*ytmp)
    do i=1,ysize(1)
      x=(i-1+ystart(1)-1)*dx
      mysx=sin(twopi*x)
      temp_top(i,j,k)=-(ytmp-2.)*fluxratio_top
    enddo
    enddo
    enddo
    !
    ! Compute current xz diffusion
    call xzdiff_temp_solide(temp_bot,bot_dm1,ny_sol_bot)
    call xzdiff_temp_solide(temp_top,top_dm1,ny_sol_top)
    ! Compute current y diffusion
    call ydiff_temp_solide(ydiff_bot,temp_bot,d2pn_bot,ny_sol_bot)
    call ydiff_temp_solide(ydiff_top,temp_top,d2pn_top,ny_sol_top)
    ! Diffusion non explicite
    bot_dm1=(bot_dm1+sol_imp_var*ydiff_bot)/repr_sol_bot
    top_dm1=(top_dm1+sol_imp_var*ydiff_top)/repr_sol_top
    if (nrank.eq.0) print *,'    - Solid temperature initialized'
    ! Compute min/max & check convergence
    if (nrank.eq.0) print *,'Solide bas : '
    call test_sol_min_max(temp_bot,yp_bot,ny_sol_bot,fluxratio_bot,repr_sol_bot)
    if (nrank.eq.0) print *,'Solide haut : '
    call test_sol_min_max(temp_top,yp_top,ny_sol_top,fluxratio_top,repr_sol_top)
    !
  else
    !
    call solide_restart(.false.)
    !
    if (nrank.eq.0) print *,'    - Solid temperature restart ok'
    !
  endif
  !
  if (nrank.eq.0) print *,'  Conjugate heat transfer initialized'
  !
end subroutine conjugate_ht_init

subroutine cheb_nodes(y,n,a,b)
!
! Store chebyshev nodes in y
!
  use param, only : pi

  implicit none

  integer, intent(in) :: n
  real(mytype), intent(in) :: a, b
  real(mytype), dimension(n+1+2), intent(out) :: y

  ! Local variables
  integer :: i

  y(1)=a
  do i=2,n+2
    y(i)=(a+b)/2.d0+(b-a)*cos(((2*(n-(i-2))+1)*1.d0/(2*n+2))*pi)/2.d0
  enddo
  y(n+1+2)=b

end subroutine cheb_nodes

subroutine cheb_der2(der2,y,n)
!
! Compute second derivative of (Chebyshev) polynoms
! At y(i) nodes
! Line for a node
! Col. for a polynom
!
  use decomp_2d, only : nproc, nrank, real_type
  use MPI

  implicit none

  integer, intent(in) :: n
  real(mytype), dimension(n+1+2), intent(in) :: y
  real(mytype), dimension(n+1,n+1), intent(out) :: der2

  ! Local variables
  integer :: i, j, k, l, code
  real(mytype), dimension(n+1) :: d2li, mysum, myprod

  ! Loop on polynoms
  if (nproc.ge.n+1) then
    !
    ! Procs 0 to n compute a column
    !
    if (nrank.le.n) then
      i=nrank+1
      d2li=0.
      do l=1,n+1
      if (l.ne.i) then
        mysum=0.
        do k=1,n+1
        if ((k.ne.i).and.(k.ne.l)) then
          myprod=1.
          do j=1,n+1
          if ((j.ne.i).and.(j.ne.k).and.(j.ne.l)) then
            myprod=myprod*(y(2:n+2)-y(j+1))/(y(i+1)-y(j+1))
          endif
          enddo
          mysum=mysum+myprod/(y(i+1)-y(k+1))
        endif
        enddo
        d2li=d2li+mysum/(y(i+1)-y(l+1))
      endif
      enddo
      der2(:,i)=d2li
    endif
    !
    ! then information is broadcasted to all procs
    !
    do i=1,n+1
      call MPI_BCAST(der2(:,i),n+1,real_type,i-1,MPI_COMM_WORLD,code)
    enddo
    !
  else
    !
    do i=1,n+1
      d2li=0.
      do l=1,n+1
      if (l.ne.i) then
        mysum=0.
        do k=1,n+1
        if ((k.ne.i).and.(k.ne.l)) then
          myprod=1.
          do j=1,n+1
          if ((j.ne.i).and.(j.ne.k).and.(j.ne.l)) then
            myprod=myprod*(y(2:n+2)-y(j+1))/(y(i+1)-y(j+1))
          endif
          enddo
          mysum=mysum+myprod/(y(i+1)-y(k+1))
        endif
        enddo
        d2li=d2li+mysum/(y(i+1)-y(l+1))
      endif
      enddo
      der2(:,i)=d2li
    enddo
    !
  endif

end subroutine cheb_der2

subroutine cheb_flux(flux,y,n)
!
! Compute First derivative of chebyshev polynoms
! At boundary nodes
! First line for derivative at first node
! Second line for derivative at last node
!
  implicit none

  integer, intent(in) :: n
  real(mytype), dimension(n+1+2), intent(in) :: y
  real(mytype), dimension(2,n+1), intent(out) :: flux

  ! Local variables
  integer :: i, j, k
  real(mytype), dimension(2) :: d1li, myprod, tmpy

  tmpy(1)=y(1)
  tmpy(2)=y(n+1+2)
  
  do i=1,n+1
    d1li=0.
    do k=1,n+1
    if (k.ne.i) then
      myprod=1.
      do j=1,n+1
      if ((j.ne.i).and.(j.ne.k)) then
        myprod=myprod*(tmpy-y(j+1))/(y(i+1)-y(j+1))
      endif
      enddo
      d1li=d1li+myprod/(y(i+1)-y(k+1))
    endif
    enddo
    flux(:,i)=d1li
  enddo

end subroutine cheb_flux

subroutine cheb_edge(val,y,n)
!
! Compute polynom value at boundary nodes
! First line for value at first node
! Second line for value at last node
!
  implicit none

  integer, intent(in) :: n
  real(mytype), dimension(n+1+2), intent(in) :: y
  real(mytype), dimension(2,n+1), intent(out) :: val

  ! Local variables
  integer :: i, j
  real(mytype), dimension(2) :: li, tmpy

  tmpy(1)=y(1)
  tmpy(2)=y(n+1+2)

  do i=1,n+1
    li=1.
    do j=1,n+1
    if (j.ne.i) then
      li=li*(tmpy-y(j+1))/(y(i+1)-y(j+1))
    endif
    enddo
    val(:,i)=li
  enddo

end subroutine cheb_edge

subroutine update_temp_solide()
!
! Compute temperature in the solid
! with the boundary flux imposed from the fluid
!
  use decomp_2d, only : ysize, nrank
  use param, only : dt, twopi, itime, ilit, t

  implicit none

  ! Local variables
  integer :: i, j, k
  real(mytype) :: a, b, c
  real(mytype), dimension(ysize(1),ny_sol_bot+1+2,ysize(3)) :: ydiff_bot
  real(mytype), dimension(ysize(1),ny_sol_top+1+2,ysize(3)) :: ydiff_top

  ! Switch from Euler1 to AB3
  if ((ilit.eq.0).and.(itime.eq.1)) then
    a=   dt
    b=   0.
    c=   0.
  elseif ((ilit.eq.0).and.(itime.eq.2)) then
    a=   1.5*dt
    b= - 0.5*dt
    c=   0.
  else
    a=   (23./12.)*dt
    b= - (16./12.)*dt
    c=   ( 5./12.)*dt
  endif

  ! Full RHS in temporary array ydiff_***
  ydiff_bot = temp_bot &
            + a*bot_dm1 &
            + b*bot_dm2 &
            + c*bot_dm3
  ydiff_top = temp_top &
            + a*top_dm1 &
            + b*top_dm2 &
            + c*top_dm3

  ! Boundary coundition : flux imposed with fluid and temperature imposed with exterior
  ! Derivate in solid is derivate in fluid * fluxratio
  ydiff_bot(:,1,:)=1.!-fluxratio_bot
  ydiff_bot(:,ny_sol_bot+1+2,:)=cl_bot(:,1,:)*fluxratio_bot
  ydiff_top(:,1,:)=cl_top(:,1,:)*fluxratio_top
  ydiff_top(:,ny_sol_top+1+2,:)=-1.!-fluxratio_top

  ! Use matrix inversion to update temperature
  ! New temperature in temp_bot and temp_top
  call my_new_solid_temp(temp_bot,ydiff_bot,invdiffy_bot,&
           d2pn_bot,val_bot,ny_sol_bot,flux_bot(1,:),flux_bot(2,:),repr_sol_bot)
  call my_new_solid_temp(temp_top,ydiff_top,invdiffy_top,&
           d2pn_top,val_top,ny_sol_top,flux_top(1,:),flux_top(2,:),repr_sol_top)

  ! Compute min/max & check convergence if required
  if (nrank.eq.0) print *,'Solide bas : '
  call test_sol_min_max(temp_bot,yp_bot,ny_sol_bot,fluxratio_bot,repr_sol_bot)
  if (nrank.eq.0) print *,'Solide haut : '
  call test_sol_min_max(temp_top,yp_top,ny_sol_top,fluxratio_top,repr_sol_top)

  ! Update history
  bot_dm3=bot_dm2
  top_dm3=top_dm2
  bot_dm2=bot_dm1
  top_dm2=top_dm1

  ! Compute current xz diffusion
  ! (bot_dm1 and top_dm1 are updated)
  call xzdiff_temp_solide(temp_bot,bot_dm1,ny_sol_bot)
  call xzdiff_temp_solide(temp_top,top_dm1,ny_sol_top)
  ! Compute current y diffusion
  ! (stored in temporary array ydiff_***)
  call ydiff_temp_solide(ydiff_bot,temp_bot,d2pn_bot,ny_sol_bot)
  call ydiff_temp_solide(ydiff_top,temp_top,d2pn_top,ny_sol_top)
  ! Diffusion non explicite
  bot_dm1=(bot_dm1+sol_imp_var*ydiff_bot)/repr_sol_bot
  top_dm1=(top_dm1+sol_imp_var*ydiff_top)/repr_sol_top

end subroutine update_temp_solide

subroutine ydiff_temp_solide(diff,temp,d2pn,n)
!
! Compute diffusion in y-direction
! Only at chebyshev nodes
!
  use decomp_2d, only : ysize

  implicit none

  integer, intent(in) :: n
  real(mytype), dimension(n+1,n+1), intent(in) :: d2pn
  real(mytype), dimension(ysize(1),n+1+2,ysize(3)), intent(in) :: temp
  real(mytype), dimension(ysize(1),n+1+2,ysize(3)), intent(out) :: diff

  ! Local variables
  integer :: j,l

  diff=0.
  do j=2,n+2
    do l=1,n+1
      diff(:,j,:)=diff(:,j,:)+d2pn(j-1,l)*temp(:,l+1,:)
    enddo
  enddo

end subroutine ydiff_temp_solide

integer function my_ydiff_mat_inv(matin,n,matout,cl_b,cl_t,repr)
!
! Compute the inverse of y diffusion
! With cl_* boundary condition
! Put it in matout
!
  use param, only : dt

  implicit none

  integer, intent(in) :: n
  real(mytype), intent(in) :: repr
  real(mytype), dimension(1,n+1), intent(in) :: cl_b, cl_t
  real(mytype), dimension(n+1,n+1), intent(in) :: matin
  real(mytype), dimension(n-1,n-1), intent(out) :: matout

  ! Local variables
  integer :: info
  integer :: i,j,k,l
  integer, dimension(n-1) :: ipiv
  real(mytype), dimension(2,2) :: math
  real(mytype), dimension(2,n-1) :: matg
  real(mytype), dimension(2,n+1) :: mycl

  real(mytype), dimension(n+1,n+1) :: mattmp
  real(mytype), dimension(:), allocatable :: work

  ! Lapack routines to inverse matrix
  external DGETRF ! LU DECOMP, double prec
  external SGETRF ! LU DECOMP, simple prec
  external DGETRI ! Matrix inversion, double prec
  external SGETRI ! Matrix inversion, simple prec

  my_ydiff_mat_inv=0
  matout=0.
  mattmp=0.

  mycl(1,:)=cl_b(1,:)
  mycl(2,:)=cl_t(1,:)

  ! Compute math (inverse of a 2x2 matrix)
  math(1,1)=mycl(2,n+1)
  math(2,2)=mycl(1,1)
  math(1,2)=-mycl(1,n+1)
  math(2,1)=-mycl(2,1)
  math=math/(mycl(1,1)*mycl(2,n+1)-mycl(1,n+1)*mycl(2,1))

  ! Compute matg
  matg=0.
  do i=1,2
  do j=1,n-1
  do l=1,2
    matg(i,j) = matg(i,j) - math(i,l)*mycl(l,j+1)
  enddo
  enddo
  enddo

  ! Compute matout (matrix to inverse)
  mattmp=0.
  do i=1,n+1
    mattmp(i,i)=1.
  enddo
  ! Diffusion is not explicit
  do i=1,n+1
  do j=1,n+1
    mattmp(i,j)=mattmp(i,j)-(1.-sol_imp_var)*dt*matin(i,j)/repr
  enddo
  enddo
  do i=1,n-1
  do j=1,n-1
    matout(i,j) = mattmp(i+1,j+1) &
                + mattmp(i+1,1)*matg(1,j) &
                + mattmp(i+1,n+1)*matg(2,j)
  enddo
  enddo

  ! Call lapack to inverse the matrix
  if (allocated(work)) deallocate(work)
  allocate(work(1))
#ifdef DOUBLE_PREC
  call DGETRF( n-1, n-1, matout, n-1, ipiv, info )
  CALL DGETRI( n-1, matout, n-1, ipiv, work, -1, info)
#else
  call SGETRF( n-1, n-1, matout, n-1, ipiv, info )
  CALL SGETRI( n-1, matout, n-1, ipiv, work, -1, info)
#endif
  i=work(1)
  info=1
  do while (info.ne.0)
    j=i
    if (allocated(work)) deallocate(work)
    allocate(work(j),stat=info)
    i=j/2
  end do
#ifdef DOUBLE_PREC
  CALL DGETRI( n-1, matout, n-1, ipiv, work, j, info)
#else
  CALL SGETRI( n-1, matout, n-1, ipiv, work, j, info)
#endif
  if (allocated(work)) deallocate(work)

  my_ydiff_mat_inv = info

  ! Terminus
  return

end function my_ydiff_mat_inv

subroutine my_new_solid_temp(temp,rhs,invdiffy,d2pn,val,n,cl_b,cl_t,repr)
!
! Update temperature temp :
! temp = invdiffy * rhs
!
  use decomp_2d, only : ysize
  use param, only : dt

  implicit none

  integer, intent(in) :: n
  real(mytype), intent(in) :: repr
  real(mytype), dimension(2,n+1), intent(in) :: val
  real(mytype), dimension(1,n+1), intent(in) :: cl_b, cl_t
  real(mytype), dimension(n+1,n+1), intent(in) :: d2pn
  real(mytype), dimension(n-1,n-1), intent(in) :: invdiffy
  real(mytype), dimension(ysize(1),n+1+2,ysize(3)), intent(in) :: rhs
  real(mytype), dimension(ysize(1),n+1+2,ysize(3)), intent(out) :: temp

  ! Local variables
  integer :: i,j,k,l
  real(mytype), dimension(2,2) :: math
  real(mytype), dimension(2,n-1) :: matg
  real(mytype), dimension(2,n+1) :: mycl

  temp=0.

  mycl(1,:)=cl_b(1,:)
  mycl(2,:)=cl_t(1,:)

  ! Compute math (inverse of a 2x2 matrix)
  math(1,1)=mycl(2,n+1)
  math(2,2)=mycl(1,1)
  math(1,2)=-mycl(1,n+1)
  math(2,1)=-mycl(2,1)
  math=math/(mycl(1,1)*mycl(2,n+1)-mycl(1,n+1)*mycl(2,1))
  ! Compute matg
  matg=0.
  do i=1,2
  do j=1,n-1
  do l=1,2
    matg(i,j) = matg(i,j) - math(i,l)*mycl(l,j+1)
  enddo
  enddo
  enddo

  ! Compute temp inside (j=3:n+1)
  do k=1,ysize(3)
  do j=3,n+1
  do i=1,ysize(1)
    do l=1,n-1
      temp(i,j,k) = temp(i,j,k) + invdiffy(j-2,l)*( rhs(i,l+2,k) &
            + ((1.-sol_imp_var)/repr)*dt*d2pn(l+1,  1)* &
                  ( math(1,1)*rhs(i,1,k)+math(1,2)*rhs(i,n+3,k) ) &
            + ((1.-sol_imp_var)/repr)*dt*d2pn(l+1,n+1)* &
                  ( math(2,1)*rhs(i,1,k)+math(2,2)*rhs(i,n+3,k) ) )
    enddo
  enddo
  enddo
  enddo

  ! Compute temp at boundary chebyshev nodes
  ! Using boundary condition
  ! j=2 & j=n+2
  !
  ! Non-homogeneous boundary-condition from fluid
  do k=1,ysize(3)
  do i=1,ysize(1)
    temp(i,2,k)   = math(1,1)*rhs(i,1,k) + math(1,2)*rhs(i,n+3,k)
    temp(i,n+2,k) = math(2,1)*rhs(i,1,k) + math(2,2)*rhs(i,n+3,k)
  enddo
  enddo
  !
  ! Homogeneous boundary condition from solid
  do k=1,ysize(3)
  do i=1,ysize(1)
    do j=1,n-1
      temp(i,2,k)   = temp(i,2,k) &
                    + matg(1,j)*temp(i,j+2,k)
      temp(i,n+2,k) = temp(i,n+2,k) & 
                    + matg(2,j)*temp(i,j+2,k)
    enddo
  enddo
  enddo

  ! Compute temp at boundary nodes
  ! Using chebyshev coefficients
  do k=1,ysize(3)
  do j=2,n+2
  do i=1,ysize(1)
    temp(i,1,k)   = temp(i,1,k) &
                  + temp(i,j,k)*val(1,j-1)
    temp(i,n+3,k) = temp(i,n+3,k) &
                  + temp(i,j,k)*val(2,j-1)
  enddo
  enddo
  enddo

end subroutine my_new_solid_temp

subroutine xzdiff_temp_solide(temp,dm1,n)
!
! Compute diffusion of temp
! In x and z direction
! and put it in dm1
!
  use decomp_2d, only : ysize, nproc, neighbour, nrank, ystart, yend, real_type
  use param, only : dx2, dz2, itime, ifirst
  use MPI

  implicit none

  integer, intent(in) :: n
  real(mytype), dimension(ysize(1),n+1+2,ysize(3)), intent(in) :: temp
  real(mytype), dimension(ysize(1),n+1+2,ysize(3)), intent(out) :: dm1

  ! Local variables
  integer status(MPI_STATUS_SIZE)
  integer :: i,j,k
  integer :: req_0s, req_0r, req_1s, req_1r
  integer, save :: x0, x1, z0, z1
  real(mytype) :: a, b, c
  real(mytype) :: ax, bx, cx
  real(mytype) :: az, bz, cz
  real(mytype), dimension(3,n+1+2,ysize(3)) :: tmpxx, tmpxx2
  real(mytype), dimension(-2:0,n+1+2,ysize(3)) :: xmoins
  real(mytype), dimension(ysize(1)+1:ysize(1)+3,n+1+2,ysize(3)) :: xplus
  real(mytype), dimension(ysize(1),n+1+2,3) :: tmpzz, tmpzz2
  real(mytype), dimension(ysize(1),n+1+2,-2:0) :: zmoins
  real(mytype), dimension(ysize(1),n+1+2,ysize(3)+1:ysize(3)+3) :: zplus

  ! Finite difference scheme for x & z diffusion
  a=6./4.
  b=-3./5.
  c=2./20.
  ax=a/dx2
  bx=b/(4.*dx2)
  cx=c/(9.*dx2)
  az=a/dz2
  bz=b/(4.*dz2)
  cz=c/(9.*dz2)

  ! Init diffusion (dm1) and boundary values
  dm1=0.
  xmoins=0.
  xplus=0.
  zmoins=0.
  zplus=0.

  if (itime.eq.ifirst) then ! Prepare communication for periodic boundary condition

    if (neighbour(2,1).ne.neighbour(2,2)) then ! x-direction

      ! X, gauche vers droite
      if (neighbour(2,2).eq.MPI_PROC_NULL) then ! first
        !send
        call MPI_SEND(nrank, 1, mpi_integer, &
                      neighbour(2,1), nrank, MPI_COMM_WORLD, i)
        call MPI_RECV(x0, 1, mpi_integer, &
                      MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, status, i)
      elseif (neighbour(2,1).ne.MPI_PROC_NULL) then ! middle
        call MPI_RECV(x0, 1, mpi_integer, &
                      neighbour(2,2), neighbour(2,2), MPI_COMM_WORLD, status, i)
        call MPI_SEND(x0, 1, mpi_integer, &
                      neighbour(2,1), nrank, MPI_COMM_WORLD, i)
        x0=MPI_PROC_NULL
      else ! last
        call MPI_RECV(x0, 1, mpi_integer, &
                      neighbour(2,2), neighbour(2,2), MPI_COMM_WORLD, status, i)
        call MPI_SEND(nrank, 1, mpi_integer, &
                      x0, 0, MPI_COMM_WORLD, i)
      endif
      call MPI_BARRIER(mpi_comm_world,i)
      ! X, droite vers gauche
      if (neighbour(2,1).eq.MPI_PROC_NULL) then ! last
        !send
        call MPI_SEND(nrank, 1, mpi_integer, &
                      neighbour(2,2), nrank, MPI_COMM_WORLD, i)
        !recv
        call MPI_RECV(x1, 1, mpi_integer, &
                      MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, status, i)
      elseif (neighbour(2,2).ne.MPI_PROC_NULL) then ! middle
        !recv
        call MPI_RECV(x1, 1, mpi_integer, &
                      neighbour(2,1), neighbour(2,1), MPI_COMM_WORLD, status, i)
        !send
        call MPI_SEND(x1, 1, mpi_integer, &
                      neighbour(2,2), nrank, MPI_COMM_WORLD, i)
        x1=MPI_PROC_NULL
      else ! first
        !recv
        call MPI_RECV(x1, 1, mpi_integer, &
                      neighbour(2,1), neighbour(2,1), MPI_COMM_WORLD, status, i)
        !send to last
        call MPI_SEND(nrank, 1, mpi_integer, &
                      x1, 1, MPI_COMM_WORLD, i)
      endif
      call MPI_BARRIER(mpi_comm_world,i)
    endif

    if (neighbour(2,5).ne.neighbour(2,6)) then ! z-direction

      ! Z, bas vers haut
      if (neighbour(2,6).eq.MPI_PROC_NULL) then ! first
        !send
        call MPI_SEND(nrank, 1, mpi_integer, &
                      neighbour(2,5), nrank, MPI_COMM_WORLD, i)
        !recv
        call MPI_RECV(z0, 1, mpi_integer, &
                      MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, status, i)
      elseif (neighbour(2,5).ne.MPI_PROC_NULL) then ! middle
        !recv
        call MPI_RECV(z0, 1, mpi_integer, &
                      neighbour(2,6), neighbour(2,6), MPI_COMM_WORLD, status, i)
        !send
        call MPI_SEND(z0, 1, mpi_integer, &
                      neighbour(2,5), nrank, MPI_COMM_WORLD, i)
        z0=MPI_PROC_NULL
      else ! last
        !recv
        call MPI_RECV(z0, 1, mpi_integer, &
                      neighbour(2,6), neighbour(2,6), MPI_COMM_WORLD, status, i)
        !send to first
        call MPI_SEND(nrank, 1, mpi_integer, &
                      z0, 2, MPI_COMM_WORLD, i)
      endif
      call MPI_BARRIER(mpi_comm_world,i)
      ! Z, haut vers bas
      if (neighbour(2,5).eq.MPI_PROC_NULL) then ! last
        !send
        call MPI_SEND(nrank, 1, mpi_integer, &
                      neighbour(2,6), nrank, MPI_COMM_WORLD, i)
        !recv
        call MPI_RECV(z1, 1, mpi_integer, &
                      MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, status, i)
      elseif (neighbour(2,6).ne.MPI_PROC_NULL) then ! middle
        !recv
        call MPI_RECV(z1, 1, mpi_integer, &
                      neighbour(2,5), neighbour(2,5), MPI_COMM_WORLD, status, i)
        !send
        call MPI_SEND(z1, 1, mpi_integer, &
                      neighbour(2,6), nrank, MPI_COMM_WORLD, i)
        z1=MPI_PROC_NULL
      else ! first
        !recv
        call MPI_RECV(z1, 1, mpi_integer, &
                      neighbour(2,5), neighbour(2,5), MPI_COMM_WORLD, status, i)
        !send to last
        call MPI_SEND(nrank, 1, mpi_integer, &
                      z1, 3, MPI_COMM_WORLD, i)
      endif
      call MPI_BARRIER(mpi_comm_world,i)
    endif

  endif

  if (neighbour(2,1).eq.neighbour(2,2)) then ! No neighbour in x-direction
    xmoins(-2,:,:)=temp(ysize(1)-2,:,:)
    xmoins(-1,:,:)=temp(ysize(1)-1,:,:)
    xmoins( 0,:,:)=temp(ysize(1)  ,:,:)
    xplus(ysize(1)+1,:,:)=temp(1,:,:)
    xplus(ysize(1)+2,:,:)=temp(2,:,:)
    xplus(ysize(1)+3,:,:)=temp(3,:,:)
  else

    ! Init
    tmpxx(:,:,:)=temp(ysize(1)-2:ysize(1),:,:)

    ! X, gauche vers droite
    if (neighbour(2,2).eq.MPI_PROC_NULL) then ! first
      !recv
      call MPI_IRECV(tmpxx2, 3*(n+1+2)*ysize(3), real_type, &
                    x0, x0, MPI_COMM_WORLD, req_0r, i)
      !send
      call MPI_ISEND(tmpxx, 3*(n+1+2)*ysize(3), real_type, &
                    neighbour(2,1), nrank, MPI_COMM_WORLD, req_0s, i)
    elseif (neighbour(2,1).ne.MPI_PROC_NULL) then ! middle
      !recv
      call MPI_IRECV(tmpxx2, 3*(n+1+2)*ysize(3), real_type, &
                    neighbour(2,2), neighbour(2,2), MPI_COMM_WORLD, req_0r, i)
      !send
      call MPI_ISEND(tmpxx, 3*(n+1+2)*ysize(3), real_type, &
                    neighbour(2,1), nrank, MPI_COMM_WORLD, req_0s, i)
    else ! last
      !recv
      call MPI_IRECV(tmpxx2, 3*(n+1+2)*ysize(3), real_type, &
                    neighbour(2,2), neighbour(2,2), MPI_COMM_WORLD, req_0r, i)
      !send to first
      call MPI_ISEND(tmpxx, 3*(n+1+2)*ysize(3), real_type, &
                    x0, nrank, MPI_COMM_WORLD, req_0s, i)
    endif
    call MPI_Wait(req_0s,status,i)
    call MPI_Wait(req_0r,status,i)

    xmoins=tmpxx2
    tmpxx(:,:,:)=temp(1:3,:,:)

    ! X, droite vers gauche
    if (neighbour(2,1).eq.MPI_PROC_NULL) then ! last
      !recv
      call MPI_IRECV(tmpxx2, 3*(n+1+2)*ysize(3), real_type, &
                    x1, x1, MPI_COMM_WORLD, req_1r, i)
      !send
      call MPI_ISEND(tmpxx, 3*(n+1+2)*ysize(3), real_type, &
                    neighbour(2,2), nrank, MPI_COMM_WORLD, req_1s, i)
    elseif (neighbour(2,2).ne.MPI_PROC_NULL) then ! middle
      !recv
      call MPI_IRECV(tmpxx2, 3*(n+1+2)*ysize(3), real_type, &
                    neighbour(2,1), neighbour(2,1), MPI_COMM_WORLD, req_1r, i)
      !send
      call MPI_ISEND(tmpxx, 3*(n+1+2)*ysize(3), real_type, &
                    neighbour(2,2), nrank, MPI_COMM_WORLD, req_1s, i)
    else ! first
      !recv
      call MPI_IRECV(tmpxx2, 3*(n+1+2)*ysize(3), real_type, &
                    neighbour(2,1), neighbour(2,1), MPI_COMM_WORLD, req_1r, i)
      !send to last
      call MPI_ISEND(tmpxx, 3*(n+1+2)*ysize(3), real_type, &
                    x1, nrank, MPI_COMM_WORLD, req_1s, i)
    endif

    call MPI_Wait(req_1s,status,i)
    call MPI_Wait(req_1r,status,i)
    xplus=tmpxx2

  endif

  if (neighbour(2,5).eq.neighbour(2,6)) then ! No neighbour in z-direction
    zmoins(:,:,-2)=temp(:,:,ysize(3)-2)
    zmoins(:,:,-1)=temp(:,:,ysize(3)-1)
    zmoins(:,:, 0)=temp(:,:,ysize(3)  )
    zplus(:,:,ysize(3)+1)=temp(:,:,1)
    zplus(:,:,ysize(3)+2)=temp(:,:,2)
    zplus(:,:,ysize(3)+3)=temp(:,:,3)
  else

    ! Init
    tmpzz(:,:,:)=temp(:,:,ysize(3)-2:ysize(3))

    ! Z, bas vers haut
    if (neighbour(2,6).eq.MPI_PROC_NULL) then ! first
      !recv
      call MPI_IRECV(tmpzz2, ysize(1)*(n+1+2)*3, real_type, &
                    z0, z0, MPI_COMM_WORLD, req_0r, i)
      !send
      call MPI_ISEND(tmpzz, ysize(1)*(n+1+2)*3, real_type, &
                    neighbour(2,5), nrank, MPI_COMM_WORLD, req_0s, i)
    elseif (neighbour(2,5).ne.MPI_PROC_NULL) then ! middle
      !recv
      call MPI_IRECV(tmpzz2, ysize(1)*(n+1+2)*3, real_type, &
                    neighbour(2,6), neighbour(2,6), MPI_COMM_WORLD, req_0r, i)
      !send
      call MPI_ISEND(tmpzz, ysize(1)*(n+1+2)*3, real_type, &
                    neighbour(2,5), nrank, MPI_COMM_WORLD, req_0s, i)
    else ! last
      !recv
      call MPI_IRECV(tmpzz2, ysize(1)*(n+1+2)*3, real_type, &
                    neighbour(2,6), neighbour(2,6), MPI_COMM_WORLD, req_0r, i)
      !send to first
      call MPI_ISEND(tmpzz, ysize(1)*(n+1+2)*3, real_type, &
                    z0, nrank, MPI_COMM_WORLD, req_0s, i)
    endif

    call MPI_Wait(req_0r,status,i)
    call MPI_Wait(req_0s,status,i)

    zmoins=tmpzz2
    tmpzz(:,:,:)=temp(:,:,1:3)

    ! Z, haut vers bas
    if (neighbour(2,5).eq.MPI_PROC_NULL) then ! last
      !recv
      call MPI_IRECV(tmpzz2, ysize(1)*(n+1+2)*3, real_type, &
                    z1, z1, MPI_COMM_WORLD, req_1r, i)
      !send
      call MPI_ISEND(tmpzz, ysize(1)*(n+1+2)*3, real_type, &
                    neighbour(2,6), nrank, MPI_COMM_WORLD, req_1s, i)
    elseif (neighbour(2,6).ne.MPI_PROC_NULL) then ! middle
      !recv
      call MPI_IRECV(tmpzz2, ysize(1)*(n+1+2)*3, real_type, &
                    neighbour(2,5), neighbour(2,5), MPI_COMM_WORLD, req_1r, i)
      !send
      call MPI_ISEND(tmpzz, ysize(1)*(n+1+2)*3, real_type, &
                    neighbour(2,6), nrank, MPI_COMM_WORLD, req_1s, i)
    else ! first
      !recv
      call MPI_IRECV(tmpzz2, ysize(1)*(n+1+2)*3, real_type, &
                    neighbour(2,5), neighbour(2,5), MPI_COMM_WORLD, req_1r, i)
      !send to last
      call MPI_ISEND(tmpzz, ysize(1)*(n+1+2)*3, real_type, &
                    z1, nrank, MPI_COMM_WORLD, req_1s, i)
    endif

    call MPI_Wait(req_1r,status,i)
    call MPI_Wait(req_1s,status,i)
    zplus=tmpzz2

  endif

  ! xdiff
  do k=1,ysize(3)
  do j=1,n+1+2
    i=1
    dm1(i,j,k) = ax*(xmoins(i-1,j,k)+temp(i+1,j,k)-2.*temp(i,j,k)) &
               + bx*(xmoins(i-2,j,k)+temp(i+2,j,k)-2.*temp(i,j,k)) &
               + cx*(xmoins(i-3,j,k)+temp(i+3,j,k)-2.*temp(i,j,k))
    i=2
    dm1(i,j,k) = ax*(temp(i-1  ,j,k)+temp(i+1,j,k)-2.*temp(i,j,k)) &
               + bx*(xmoins(i-2,j,k)+temp(i+2,j,k)-2.*temp(i,j,k)) &
               + cx*(xmoins(i-3,j,k)+temp(i+3,j,k)-2.*temp(i,j,k))
    i=3
    dm1(i,j,k) = ax*(temp(i-1  ,j,k)+temp(i+1,j,k)-2.*temp(i,j,k)) &
               + bx*(temp(i-2  ,j,k)+temp(i+2,j,k)-2.*temp(i,j,k)) &
               + cx*(xmoins(i-3,j,k)+temp(i+3,j,k)-2.*temp(i,j,k))
    do i=4,ysize(1)-3
      dm1(i,j,k) = ax*(temp(i-1,j,k)+temp(i+1,j,k)-2.*temp(i,j,k)) &
                 + bx*(temp(i-2,j,k)+temp(i+2,j,k)-2.*temp(i,j,k)) &
                 + cx*(temp(i-3,j,k)+temp(i+3,j,k)-2.*temp(i,j,k))
    enddo
    i=ysize(1)-2
    dm1(i,j,k) = ax*(temp(i-1,j,k)+temp(i+1 ,j,k)-2.*temp(i,j,k)) &
               + bx*(temp(i-2,j,k)+temp(i+2 ,j,k)-2.*temp(i,j,k)) &
               + cx*(temp(i-3,j,k)+xplus(i+3,j,k)-2.*temp(i,j,k))
    i=ysize(1)-1
    dm1(i,j,k) = ax*(temp(i-1,j,k)+temp(i+1 ,j,k)-2.*temp(i,j,k)) &
               + bx*(temp(i-2,j,k)+xplus(i+2,j,k)-2.*temp(i,j,k)) &
               + cx*(temp(i-3,j,k)+xplus(i+3,j,k)-2.*temp(i,j,k))
    i=ysize(1)
    dm1(i,j,k) = ax*(temp(i-1,j,k)+xplus(i+1,j,k)-2.*temp(i,j,k)) &
               + bx*(temp(i-2,j,k)+xplus(i+2,j,k)-2.*temp(i,j,k)) &
               + cx*(temp(i-3,j,k)+xplus(i+3,j,k)-2.*temp(i,j,k))
  enddo
  enddo

  ! + zdiff
  k=1
  do j=1,n+1+2
  do i=1,ysize(1)
    dm1(i,j,k) = dm1(i,j,k) &
               + az*(zmoins(i,j,k-1)+temp(i,j,k+1)-2.*temp(i,j,k)) &
               + bz*(zmoins(i,j,k-2)+temp(i,j,k+2)-2.*temp(i,j,k)) &
               + cz*(zmoins(i,j,k-3)+temp(i,j,k+3)-2.*temp(i,j,k))
  enddo
  enddo
  k=2
  do j=1,n+1+2
  do i=1,ysize(1)
    dm1(i,j,k) = dm1(i,j,k) &
               + az*(temp(i,j,k-1  )+temp(i,j,k+1)-2.*temp(i,j,k)) &
               + bz*(zmoins(i,j,k-2)+temp(i,j,k+2)-2.*temp(i,j,k)) &
               + cz*(zmoins(i,j,k-3)+temp(i,j,k+3)-2.*temp(i,j,k))
  enddo
  enddo
  k=3
  do j=1,n+1+2
  do i=1,ysize(1)
    dm1(i,j,k) = dm1(i,j,k) &
               + az*(temp(i,j,k-1  )+temp(i,j,k+1)-2.*temp(i,j,k)) &
               + bz*(temp(i,j,k-2  )+temp(i,j,k+2)-2.*temp(i,j,k)) &
               + cz*(zmoins(i,j,k-3)+temp(i,j,k+3)-2.*temp(i,j,k))
  enddo
  enddo
  do k=4,ysize(3)-3
  do j=1,n+1+2
  do i=1,ysize(1)
    dm1(i,j,k) = dm1(i,j,k) &
               + az*(temp(i,j,k-1)+temp(i,j,k+1)-2.*temp(i,j,k)) &
               + bz*(temp(i,j,k-2)+temp(i,j,k+2)-2.*temp(i,j,k)) &
               + cz*(temp(i,j,k-3)+temp(i,j,k+3)-2.*temp(i,j,k))
  enddo
  enddo
  enddo
  k=ysize(3)-2
  do j=1,n+1+2
  do i=1,ysize(1)
    dm1(i,j,k) = dm1(i,j,k) &
               + az*(temp(i,j,k-1)+temp(i,j,k+1 )-2.*temp(i,j,k)) &
               + bz*(temp(i,j,k-2)+temp(i,j,k+2 )-2.*temp(i,j,k)) &
               + cz*(temp(i,j,k-3)+zplus(i,j,k+3)-2.*temp(i,j,k))
  enddo
  enddo
  k=ysize(3)-1
  do j=1,n+1+2
  do i=1,ysize(1)
    dm1(i,j,k) = dm1(i,j,k) &
               + az*(temp(i,j,k-1)+temp(i,j,k+1 )-2.*temp(i,j,k)) &
               + bz*(temp(i,j,k-2)+zplus(i,j,k+2)-2.*temp(i,j,k)) &
               + cz*(temp(i,j,k-3)+zplus(i,j,k+3)-2.*temp(i,j,k))
  enddo
  enddo
  k=ysize(3)
  do j=1,n+1+2
  do i=1,ysize(1)
    dm1(i,j,k) = dm1(i,j,k) &
               + az*(temp(i,j,k-1)+zplus(i,j,k+1)-2.*temp(i,j,k)) &
               + bz*(temp(i,j,k-2)+zplus(i,j,k+2)-2.*temp(i,j,k)) &
               + cz*(temp(i,j,k-3)+zplus(i,j,k+3)-2.*temp(i,j,k))
  enddo
  enddo

end subroutine xzdiff_temp_solide

subroutine test_sol_min_max(temp,y,n,ratio,repr)
!
! Compute minimum & maximum value on both sides of the solid
! Can be used to check convergence with objective function
!
  use decomp_2d, only : ysize, nrank, real_type, ystart
  use param, only : t, twopi, dt, dx, dz
  use variables, only : nx, nz
  use MPI
!
  implicit none
!
  integer, intent(in) :: n
  real(mytype), intent(in) :: ratio, repr
  real(mytype), dimension(n+1+2), intent(in) :: y
  real(mytype), dimension(ysize(1),n+1+2,ysize(3)), intent(in) :: temp

  ! Local variables
  integer :: i,j,k,l
  integer, dimension(3) :: mymax
  real(mytype) :: x, z, mysx, mysy, mysz
  real(mytype) :: phimax, phimin, phimax1, phimin1
  real(mytype), dimension(ysize(1),n+1+2,ysize(3)) :: temp_obj

  ! Print min/max
  phimax=maxval(temp)
  phimin=minval(temp)
  call MPI_REDUCE(phimax,phimax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,i)
  call MPI_REDUCE(phimin,phimin1,1,real_type,MPI_MIN,0,MPI_COMM_WORLD,i)
  if (nrank.eq.0) print *,'Temperature min-max : ',phimin1, phimax1

  ! Check convergence on temperature
  if (.false.) then
    do k=1,ysize(3)
      z=(k-1+ystart(3)-1)*dz
      mysz=sin(twopi*z)
    do j=1,n+1+2
      mysy=sin(twopi*y(j))
    do i=1,ysize(1)
      x=(i-1+ystart(1)-1)*dx
      mysx=sin(twopi*x)
      temp_obj(i,j,k)= mysx*mysy*mysz*ratio*exp(-(twopi**2)*t*3./repr)
    enddo
    enddo
    enddo
    ! Norme L_infinie
    phimax=maxval(abs(temp-temp_obj))
    call MPI_ALLREDUCE(phimax,phimax1,1,real_type,MPI_MAX,MPI_COMM_WORLD,i)
    if (nrank.eq.0) print *,'Norme T infinie = ',phimax1
    if (phimax.eq.phimax1) then
      mymax=maxloc(abs(temp-temp_obj))
      print *,' at ',mymax+ystart,' on ',nrank
    endif
    ! Norme L2
    phimax=sum((temp-temp_obj)**2)
    call MPI_REDUCE(phimax,phimax1,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,i)
    if (nrank.eq.0) phimax1=phimax1/(nx*nz*(n+1+2))
    if (nrank.eq.0) print *,'Norme T2 (**1/2)= ',sqrt(phimax1)

!    print *,'temperature calcul : '
!    print *,(temp(1,i,1),i=1,ny_sol_bot+3)

!    print *,'temperature delta : '
!    print *,(temp(1,i,1)-temp_obj(1,i,1),i=1,ny_sol_bot+3)
!    print *,maxloc(abs(temp-temp_obj))
  endif

end subroutine test_sol_min_max

subroutine solide_restart(will_write)
!
! Subroutine called to read or write restart file
! Use MPI-IO and 2DECOMP to read / write
!
  use decomp_2d_io, only : decomp_2d_write_var, decomp_2d_read_var, decomp_2d_write_one
  use param, only : itime, isave
  use MPI

  implicit none

  logical, intent(in) :: will_write

  ! Local variable
  integer :: fh, ierror
  integer (kind=MPI_OFFSET_KIND) :: filesize, disp
  character(len=3) :: compteur

  if (will_write) then ! write restart file
    call MPI_FILE_OPEN(MPI_COMM_WORLD, 'solide.dat', &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
         fh, ierror)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
    call decomp_2d_write_var(fh,disp,2,temp_bot,mydecomp_bot)
    call decomp_2d_write_var(fh,disp,2,bot_dm1 ,mydecomp_bot)
    call decomp_2d_write_var(fh,disp,2,bot_dm2 ,mydecomp_bot)
    call decomp_2d_write_var(fh,disp,2,bot_dm3 ,mydecomp_bot)
    call decomp_2d_write_var(fh,disp,2,temp_top,mydecomp_top)
    call decomp_2d_write_var(fh,disp,2,top_dm1 ,mydecomp_top)
    call decomp_2d_write_var(fh,disp,2,top_dm2 ,mydecomp_top)
    call decomp_2d_write_var(fh,disp,2,top_dm3 ,mydecomp_top)
    call MPI_FILE_CLOSE(fh,ierror)
    write(compteur,'(I3.3)') itime/isave
    call decomp_2d_write_one(2,temp_bot,'temp_bot'//compteur,mydecomp_bot)
    call decomp_2d_write_one(2,temp_top,'temp_top'//compteur,mydecomp_top)
  else ! read restart file
    call MPI_FILE_OPEN(MPI_COMM_WORLD, 'solide.dat', &
         MPI_MODE_RDONLY, MPI_INFO_NULL, &
         fh, ierror)
    disp = 0_MPI_OFFSET_KIND
    call decomp_2d_read_var(fh,disp,2,temp_bot,mydecomp_bot)
    call decomp_2d_read_var(fh,disp,2,bot_dm1 ,mydecomp_bot)
    call decomp_2d_read_var(fh,disp,2,bot_dm2 ,mydecomp_bot)
    call decomp_2d_read_var(fh,disp,2,bot_dm3 ,mydecomp_bot)
    call decomp_2d_read_var(fh,disp,2,temp_top,mydecomp_top)
    call decomp_2d_read_var(fh,disp,2,top_dm1 ,mydecomp_top)
    call decomp_2d_read_var(fh,disp,2,top_dm2 ,mydecomp_top)
    call decomp_2d_read_var(fh,disp,2,top_dm3 ,mydecomp_top)
    call MPI_FILE_CLOSE(fh,ierror)
  endif

end subroutine solide_restart

subroutine update_solide_stats()
!
!  Subroutine updating statistics
!
  use decomp_2d, only : nrank
  use param, only : itime

  implicit none

  ! Local variables
  integer :: coeff

  coeff=itime-beg_stat_sol
  if (coeff.le.0) then
    if (nrank.eq.0) print*,'ERROR : cannot start to compute solid statistics now, skip.'
    return
  endif

  tempm_1_bot=((coeff-1)*1./coeff)*tempm_1_bot + (1./coeff)*temp_bot
  tempm_2_bot=((coeff-1)*1./coeff)*tempm_2_bot + (1./coeff)*temp_bot**2

  tempm_1_top=((coeff-1)*1./coeff)*tempm_1_top + (1./coeff)*temp_top
  tempm_2_top=((coeff-1)*1./coeff)*tempm_2_top + (1./coeff)*temp_top**2

end subroutine update_solide_stats

subroutine solide_stats_restart(will_write)
!
!  Subroutine to read/write solid statistics
!
  use decomp_2d_io, only : decomp_2d_write_one, decomp_2d_read_one
  use MPI

  implicit none

  logical, intent(in) :: will_write

  ! Local variable
  logical :: file_exist
  integer :: fh, ierror
  integer (kind=MPI_OFFSET_KIND) :: filesize, disp

  if (will_write) then ! write restart file
    call decomp_2d_write_one(2,tempm_1_bot,'sol_bot_tempm.dat' ,mydecomp_bot)
    call decomp_2d_write_one(2,tempm_2_bot,'sol_bot_temp2m.dat',mydecomp_bot)
    call decomp_2d_write_one(2,tempm_1_top,'sol_top_tempm.dat' ,mydecomp_top)
    call decomp_2d_write_one(2,tempm_2_top,'sol_top_temp2m.dat',mydecomp_top)
  else ! read restart file
    INQUIRE(FILE='sol_bot_tempm.dat', EXIST=file_exist)
    if (file_exist) then
      call decomp_2d_read_one(2,tempm_1_bot,'sol_bot_tempm.dat' ,mydecomp_bot)
    else
      tempm_1_bot=0.
    endif
    INQUIRE(FILE='sol_bot_temp2m.dat', EXIST=file_exist)
    if (file_exist) then
      call decomp_2d_read_one(2,tempm_2_bot,'sol_bot_temp2m.dat',mydecomp_bot)
    else
      tempm_2_bot=0.
    endif
    INQUIRE(FILE='sol_top_tempm.dat', EXIST=file_exist)
    if (file_exist) then
      call decomp_2d_read_one(2,tempm_1_top,'sol_top_tempm.dat',mydecomp_top)
    else
      tempm_1_top=0.
    endif
    INQUIRE(FILE='sol_top_temp2m.dat', EXIST=file_exist)
    if (file_exist) then
      call decomp_2d_read_one(2,tempm_2_top,'sol_top_temp2m.dat',mydecomp_top)
    else
      tempm_2_top=0.
    endif
  endif

end subroutine solide_stats_restart

subroutine allocate_solide_stats()
!
!  Subroutine to allocate solid statistics
!
  use decomp_2d, only : ysize

  implicit none

  ! bottom
  allocate( tempm_1_bot(ysize(1),ny_sol_bot+1+2,ysize(3)) )
  allocate( tempm_2_bot(ysize(1),ny_sol_bot+1+2,ysize(3)) )
  tempm_1_bot=0.
  tempm_2_bot=0.

  ! top
  allocate( tempm_1_top(ysize(1),ny_sol_top+1+2,ysize(3)) )
  allocate( tempm_2_top(ysize(1),ny_sol_top+1+2,ysize(3)) )
  tempm_1_top=0.
  tempm_2_top=0.

end subroutine allocate_solide_stats

subroutine my2decomp_solide()
!
! 2DECOMP information ONLY for Y-stencil MPI IO
!
  use decomp_2d, only : ysize, ystart, get_decomp_info

  implicit none

  type(decomp_info) :: decomp_main

  call get_decomp_info(decomp_main)

  ! Define solid 2DECOMP Y-stencil information
  ! xsz(1), ysz(2), zsz(3) : total domain size
  mydecomp_bot%xsz(1)=decomp_main%xsz(1)
  mydecomp_bot%ysz(2)=ny_sol_bot+1+2
  mydecomp_bot%zsz(3)=decomp_main%zsz(3)
  mydecomp_top%xsz(1)=decomp_main%xsz(1)
  mydecomp_top%ysz(2)=ny_sol_top+1+2
  mydecomp_top%zsz(3)=decomp_main%zsz(3)
  ! ysz(1), ysz(3) : local domain size
  mydecomp_bot%ysz(1)=ysize(1)
  mydecomp_bot%ysz(3)=ysize(3)
  mydecomp_top%ysz(1)=ysize(1)
  mydecomp_top%ysz(3)=ysize(3)
  ! yst(1), yst(2), yst(3) : local domain location
  mydecomp_bot%yst(1)=ystart(1)
  mydecomp_bot%yst(2)=1
  mydecomp_bot%yst(3)=ystart(3)
  mydecomp_top%yst(1)=ystart(1)
  mydecomp_top%yst(2)=1
  mydecomp_top%yst(3)=ystart(3)

end subroutine my2decomp_solide

end module conjugate_ht
