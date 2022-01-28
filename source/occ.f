      subroutine occ(ndim,nel,nopen,ihomo,na,nb,focc)
      implicit none
      integer nel,nopen,ndim,ihomo,na,nb
      real*8 focc(ndim)
      integer i

      focc=0
c even nel      
      if(mod(nel,2).eq.0)then
      na=nel/2
      nb=nel/2
      ihomo=nel/2
      do i=1,ihomo 
         focc(i)=2.0d0
      enddo
      if(2*ihomo.ne.nel) then
         ihomo=ihomo+1
         focc(ihomo)=1.0d0
         if(nopen.eq.0)nopen=1
      endif
      if(nopen.gt.1)then
         do i=1,nopen/2
            focc(ihomo-i+1)=focc(ihomo-i+1)-1.0
            focc(ihomo+i)=focc(ihomo+i)+1.0
         enddo
         na=na+nopen/2
         nb=nb-nopen/2
      endif
c odd nel      
      else
      na=nel/2+(nopen-1)/2+1
      nb=nel/2-(nopen-1)/2
      do i=1,na             
         focc(i)=focc(i)+1.
      enddo
      do i=1,nb             
         focc(i)=focc(i)+1.
      enddo
      endif

      do i=1,ndim
         if(focc(i).gt.0.99) ihomo=i
      enddo

      end

