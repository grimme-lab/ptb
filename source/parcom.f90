module parcom
      use iso_fortran_env, only : wp => real64
      implicit none
      private :: wp
      public

      real(wp) glob_par  (20)     
      real(wp) expscal (4,10,86)

      real(wp) ener_par1 (10,86)
      real(wp) ener_par2 (10,86)

      real(wp) shell_xi  (10,86)
      real(wp) shell_cnf1(10,86)
      real(wp) shell_cnf2(10,86)
      real(wp) shell_cnf3(10,86)
      real(wp) shell_cnf4(10,86)
      real(wp) shell_resp(10,86,2)

!     real(wp),parameter :: mull_loew14 = 0.1666666667_wp ! Mulliken-Loewdin mixing factor (0=M, 0.5=L)
      real(wp),parameter :: mull_loew14 = 0.2500000000_wp ! Mulliken-Loewdin mixing factor (0=M, 0.5=L)
                                                          ! was 1/4 for quite some time

end module parcom       
