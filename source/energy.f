      subroutine energy_hf(ndim,F,h,P,eel)
      implicit none
      integer ndim
      real*8 eel
      real*8 F(ndim*(ndim+1)/2)
      real*8 h(ndim*(ndim+1)/2)
      real*8 P(ndim*(ndim+1)/2)

      integer i,j,ij

      eel=0.0d0
      ij=0
      do i=1,ndim
         do j=1,i-1
            ij = ij +1
            eel=eel+P(ij)*(F(ij)+h(ij))
         enddo
         ij = ij +1
         eel=eel+P(ij)*(F(ij)+h(ij))*0.5d0
      enddo

      end

      subroutine energy(ndim,H0,P,eel)
      implicit none
      integer ndim
      real*8 P (ndim*(ndim+1)/2)
      real*8 H0(ndim*(ndim+1)/2)
      real*8 eel

      integer i,j,ij

      eel=0.0d0
      ij=0
      do i=1,ndim
         do j=1,i-1
            ij = ij +1
            eel=eel+P(ij)*H0(ij)
         enddo
         ij = ij +1
         eel=eel+P(ij)*H0(ij)*0.5d0
      enddo
      eel = eel * 2.0d0

      end
