subroutine sscal(n,at,psh,p1,p2,scal)
   use iso_fortran_env, only : wp => real64
   use parcom
   use bascom
   implicit none 
!! ------------------------------------------------------------------------
!  Input
!! ------------------------------------------------------------------------
   integer, intent(in)    :: n                  ! number of atoms 
   integer, intent(in)    :: at(n)              ! ordinal number of atoms
   real(wp),intent(in)    ::psh(10,n)           ! shell populations       
   real(wp),intent(in)    ::p1,p2               ! parameters               
!! ------------------------------------------------------------------------
   real(wp),intent(out)   ::scal(10,n)          ! exponent scaling factors
!! ------------------------------------------------------------------------
   integer  :: i,ish
   real(wp) :: atocc(10), qa

   do i=1,n
      call shellocc_ref(at(i),atocc) ! ref. atomic pop.
      do ish=1,bas_nsh(at(i))
         qa = atocc(ish)-psh(ish,i)
         scal(ish,i) = expscal(1,ish,at(i)) + p1*qa + p2*qa**2 ! + cn(i)*ener_par1(6,at(i))*0.01
!        scal(ish,i) = expscal(1,ish,at(i)) + ener_par1(6,at(i))*qa + p2*qa**2
         if(scal(ish,i).lt.0.05) scal(ish,i)=0.05     
         if(scal(ish,i).gt.20.0) scal(ish,i)=20.0      
      enddo
   enddo

end

