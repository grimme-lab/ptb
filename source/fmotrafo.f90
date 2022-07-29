subroutine fmotrafo(ndim,homo,F,U)
   use iso_fortran_env, only : wp => real64
   use parcom
   implicit none

!! ------------------------------------------------------------------------
!  Input
!! ------------------------------------------------------------------------
   integer, intent(in)    :: ndim               ! number of SAOs       
   integer, intent(in)    :: homo               ! as the name says       
   real(wp),intent(in)    :: F(ndim*(ndim+1)/2) ! Fock matrix
   real(wp),intent(in)    :: U(ndim,ndim)       ! MOs            
!! ------------------------------------------------------------------------
!  output (to STDIO)
!! ------------------------------------------------------------------------
!  real(wp),intent(out)   :: eps(ndim)          ! corrected orbital energies
!! ------------------------------------------------------------------------
   real*4,allocatable   :: tmp(:,:), C(:,:), FF(:,:) 
   real*4               :: fmaxia,fval
   integer              :: i,j

   allocate(tmp(ndim,ndim),C(ndim,ndim),FF(ndim,ndim))

   call blowsym84(ndim,F,FF)

   C = U  ! convert to real*4 

   call ssymm('L','L',ndim,ndim,1E0,FF,ndim,C,ndim,0E0,tmp,ndim)
   call sgemm('T','N',ndim,ndim,ndim,1E0,C,ndim,tmp,ndim,0E0,FF,ndim)

!  call prmat4(6,FF,ndim,ndim,'Fmo')

!  write(*,*) 'corrected gap ',0.85*27.2113957d0*(ff(homo+1,homo+1)-ff(homo,homo))
!  do i=1,ndim
!     eps(i)= glob_par(15) + FF(i,i) + glob_par(17)*FF(i,i)
!  enddo

   fmaxia = 0
   do i=1, homo
      do j=homo+1,ndim
          fval = abs(FF(j,i)) / (0.20+ FF(j,j) - FF(i,i)) ! correlates perfectly with molecule "difficulty"
          if(fval.gt.fmaxia) then                         ! regularized to avoid overshooting for very small gap systems
             fmaxia=fval                                  ! i.e., 0.01 for alkanes, 0.02 for alkenes, > 0.03-0.04 for TM compounds
                                                          ! difficult, low gap TM compounds can have values > 0.1
          endif
      enddo
   enddo

!  fmaxia = 0
!  do i=1, homo
!     do j=homo+1,ndim
!        fmaxia = fmaxia + abs(FF(j,i)) / abs(FF(j,j) - FF(i,i))
!     enddo
!  enddo
!  fmaxia = fmaxia / dble(homo)**1.25 ! makes it indep of mol size (but is inconsistent for folded vs. lin alkanes)
!  damp = 0.55+glob_par(9)/(1.0d0+(glob_par(10)/fmaxia)**8)  ! not used anymore, damp=0.6 fixed which is simpler and not really worse

   write(*,'('' gap in "third" F          : '',2f9.3)') 27.2113957d0*(ff(homo+1,homo+1)-ff(homo,homo))*(1d0+glob_par(17))
   write(*,'('' max FMOia                 : '',2E16.8)') fmaxia

   end
