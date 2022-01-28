subroutine CalcGammaF(F,t,m,ACCY)
  implicit none
  integer m,i
  real*8 t,accy,F(0:m)
  real*8 PT5PM,t2,exp_t,term,summ
  real*8 val(1:m+1)

! TM version
  call gammav(t,m,val)
  F(0:m)=val(1:m+1)
  return

! PT5PM=float(m)+0.5

! if (t>30.0)then
!   F(0)= 0.88622692545275801365D0/sqrt(t) ! Sqrt(pi)/2
!   do i=1,m
!      F(i) = F(i-1)*(float(i)-0.5)/t
!   enddo
!   return
! else
!   t2   = 2.0*t
!   exp_t= exp(-t)
!   term = 0.5/PT5PM
!   summ = term
!   do i=1,1000
!     term = term*t/(PT5PM+float(i))
!     summ = summ+ term
!     if (term<ACCY) goto 999
!   enddo
!999 continue       
!   F(m)= exp_t*summ
!   do i=m-1,0,-1
!    F(i)=(F(i+1)*t2+exp_t)/(2.0*float(i)+1.0)
!   enddo
! endif

end

