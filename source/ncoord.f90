
subroutine ncoord_erf(nat,at,rab,kn,cn)
   use iso_fortran_env, only : wp => real64
   use com, only : rcov
   use parcom

   implicit none

   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: rab(nat*(nat+1)/2)
   real(wp),intent(in)  :: kn         
   real(wp),intent(out) :: cn(nat)

   integer  :: i, j, k
   real(wp) :: r, rcovij, tmp, arg

   cn  = 0._wp

   k = 0
   do i = 1, nat
      do j = 1, i-1
         k = k + 1
         rcovij=(4./3.)*(rcov(at(i))+rcov(at(j)))
         arg = (rab(k)-rcovij)/rcovij
         tmp = 0.5d0 * (1d0 + erf(kn*arg)) 
         cn(i) = cn(i) + tmp      
         cn(j) = cn(j) + tmp 
      enddo
      k = k + 1
   enddo

end subroutine ncoord_erf
