program modfitpar 
   
      implicit none
      integer i,ip(150,88)

      open(unit=1,file='.fitpar')
      read(1,'(150i1)') ip(1:150,87)
      read(1,'(150i1)') ip(1:150,88)
      do i=1,86
         read(1,'(150i1)') ip(1:150,i)
      enddo
      close(1)

      open(unit=1,file='.fitpar')
      write(1,'(150i1)') ip(1:150,87)
      write(1,'(150i1)') ip(1:150,88)
      do i=1,86
         if(sum(ip(1:150,i)).gt.0)then
         ip(89,i)=1
         ip(99,i)=1
         ip(127,i)=0
         ip(128,i)=0
         endif
         write(1,'(150i1)') ip(1:150,i)
      enddo
      close(1)

End

