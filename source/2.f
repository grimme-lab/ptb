 10   read(*,*) a
!     write(*,*) ((a-b)/(a+b))**2
!     write(*,*) sqrt(a*b) ! /(0.5*(a+b))
!     write(*,*) 2/(1/a+1/b)              
!     goto 10
      write(*,*) sign(sqrt(abs(a)),a)
      end
