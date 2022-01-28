!! ------------------------------------------------------------------------
!  modify electric field perturbation
!! ------------------------------------------------------------------------

subroutine modfield(n,ndim,at,z,pa,d)
   use iso_fortran_env, only : wp => real64
   use parcom
   implicit none 
!! ------------------------------------------------------------------------
!  Input
!! ------------------------------------------------------------------------
   integer, intent(in)    :: n                  ! number of atoms 
   integer, intent(in)    :: ndim               ! number of AOs       
   integer, intent(in)    :: at(n)              ! ordinal number of atoms
   real(wp),intent(in)    :: z(n)               ! nuclear charges          
   real(wp),intent(in)    ::pa(n)               ! atomic populations       
   real(wp),intent(in)    ::d(ndim*(ndim+1)/2,3)! field ints
!! ------------------------------------------------------------------------
   integer :: i,j,m,ij,ia,ib
   real(wp) :: qa,qb

   do m=1,3
   ij = 0
   do i=1,ndim
      ia = aoat(i)
      qa = z(ia)-pa(ia)
      fa = shell_xi(9,at(ia))
      do j=1,i  
         ib = aoat(j)
         qb = z(ib)-pa(ib)
         fb = shell_xi(9,at(ib))
         ij = ij + 1
         d(ij,m) = d(ij,m) * (1_wp + qa*fa + qb*fb )                 
      enddo
   enddo
   enddo

end
