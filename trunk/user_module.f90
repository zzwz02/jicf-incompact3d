#define my_mod_stats
!#define my_mod_solide

! #ifdef my_mod_solide
! #include "myconjht.f90"
! #endif

#ifdef my_mod_stats
#include "mystats.f90"
#endif

module user_specific
!
! Module reserved for user_specific functions
!

contains

subroutine module_user_init(phG,ph1,ph2,ph3,ph4)
  !
  ! Call before main time loop
  !
  use decomp_2d, only : nrank, nx_global, ny_global, nz_global, DECOMP_INFO, mytype
  use param, only : ifirst, ilast, xlx, yly, zlz, ilit

#ifdef my_mod_stats
  use user_stats
#endif
#ifdef my_mod_solide
  use conjugate_ht
  use param, only : xnu, sc, iscalar
#endif

  implicit none

  TYPE(DECOMP_INFO), intent(in) :: phG,ph1,ph2,ph3,ph4

  integer :: j
!  real(mytype) :: xx, yy, zz

bool_conjugate_ht=.false.
#ifdef my_mod_stats
  bool_user_stat=.true.
#endif
#ifdef my_mod_solide
  if (iscalar.eq.1) bool_conjugate_ht=.true.
  if (bool_conjugate_ht) then
    ny_sol_bot=64
    ly_sol_bot=1.
    repr_sol_bot=sc/xnu
    fluxratio_bot=1.
    ny_sol_top=ny_sol_bot
    ly_sol_top=ly_sol_bot
    repr_sol_top=repr_sol_bot
    fluxratio_top=fluxratio_bot
    call conjugate_ht_init()
  endif
#endif

#ifdef my_mod_stats
  if (bool_user_stat) then
    beg_stat=10000000
    call allocate_user_stats(nx_global, ny_global, nz_global,phG,ph1,ph2,ph3,ph4)
    call read_user_stats(phG,ph1,ph2,ph3,ph4)
  endif
#endif

#ifdef my_mod_solide
  if (bool_conjugate_ht) then
#ifdef my_mod_stats
    bool_sol_stats=bool_user_stat
#else
    bool_sol_stats=.false.
#endif
    if (bool_sol_stats) then
      beg_stat_sol=beg_stat
      call allocate_solide_stats()
      call solide_stats_restart(.false.)
    endif
  endif
#endif

return
end subroutine module_user_init

subroutine module_user_write(phG,ph1,ph2,ph3,ph4)
  !
  ! Call at the end of every time step
  !
  use decomp_2d, only : DECOMP_INFO
  use param, only : itime

#ifdef my_mod_stats
  use user_stats, only : bool_user_stat, beg_stat, &
                         pre_update_user_stats, plot_my_stats
#endif
#ifdef my_mod_solide
  use conjugate_ht, only : bool_conjugate_ht, bool_sol_stats, &
                           update_solide_stats, beg_stat_sol, &
                           temp_bot
#endif

  implicit none

  TYPE(DECOMP_INFO), intent(in) :: phG,ph1,ph2,ph3,ph4

#ifdef my_mod_stats
  if (bool_user_stat) then
    if (itime.gt.beg_stat) then
      call pre_update_user_stats(phG,ph1,ph2,ph3,ph4)
    endif
      call plot_my_stats()
  endif
#endif

#ifdef my_mod_solide
  if (bool_sol_stats) then
    if (itime.gt.beg_stat_sol) call update_solide_stats()
  endif
#endif

return
end subroutine module_user_write

subroutine module_user_post(phG,ph1,ph2,ph3,ph4)
  !
  ! Call at the end of the program
  !
  use decomp_2d, only : DECOMP_INFO

#ifdef my_mod_stats
  use user_stats, only : bool_user_stat, write_user_stats
#endif
#ifdef my_mod_solide
  use conjugate_ht, only : bool_conjugate_ht, solide_restart, &
                           bool_sol_stats, solide_stats_restart
#endif

  implicit none

  TYPE(DECOMP_INFO), intent(in) :: phG,ph1,ph2,ph3,ph4

#ifdef my_mod_stats
  if (bool_user_stat) call write_user_stats(phG,ph1,ph2,ph3,ph4)
#endif

#ifdef my_mod_solide
  if (bool_conjugate_ht) then
    call solide_restart(.true.)
    if (bool_sol_stats) then
      call solide_stats_restart(.true.)
    endif
  endif
#endif

return
end subroutine module_user_post

end module user_specific
