      integer valel(86)
      data valel / 
     .1,                                               2, !He 
     .1,2,                                   3,4,5,6,7,8, !Ne 
     .1,2,                                   3,4,5,6,7,8, !Ar 
     .1,2,11,     12,13,14,15,16,17,18,19,20,3,4,5,6,7,8, !Kr 
     .1,2,11,     12,13,14,15,16,17,18,19,20,3,4,5,6,7,8, !Xe 
     .1,2,11,14*0,12,13,14,15,16,17,18,19,20,3,4,5,6,7,8/
      integer i,new,n,nl,list(86)
      real*8  ,allocatable :: xyz(:,:)
      integer, allocatable :: at(:)

      call rd0('coord',n)
      allocate(at(n),xyz(3,n))
      call rd(.false.,'coord',n,xyz,at)

!     call system('rm -rf newref')
!     new=0
!     do i=1,n
!        if(valel(at(i)).lt.2.and.at(i).gt.2) new=new+1
!     enddo
!     write(*,*) n,new
!     if(new.gt.0) call system('touch newref')

      open(unit=1,file='~/.elist')
      read(1,*) nl
      read(1,*) list(1:nl)
      close(1)

      new=0
      do i=1,n
         do j=1,nl
            if(at(i).eq.list(j)) new = new + 1
         enddo
      enddo
      if(new.ne.0) then
!        call system('pwd')
         call system('pwd; rmec; egtb')
      endif

 10   continue     

      end
