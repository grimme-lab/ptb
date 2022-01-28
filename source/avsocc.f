      real*8 socc(10000,10,86),aocc(10),tmp
      integer nmol(86)

      socc = 0
      nmol = 0
      open(unit=1,file='.data')
10    read(1,*,end=99) j, i, tmp
      call getnsh(j,nsh)
      nmol(j)=nmol(j)+1
      socc(nmol(j),1,j)=tmp
      do i=2,nsh
         read(1,*,end=99) jj, ii, socc(nmol(j),i,j)
      enddo
      goto 10
99    continue
      
      do j=1,86
         if(nmol(j).eq.0) cycle
         call getnsh(j,nsh)
         aocc = 0
         do k=1,nmol(j)
!           write(*,'(i4,5f8.3)') k,socc(k,1:nsh,j)
            aocc(1:nsh)=aocc(1:nsh)+socc(k,1:nsh,j)
         enddo
         aocc = aocc /dble(nmol(j))
         write(*,*) j,nmol(j),sum(aocc)
         do i=1,nsh
            write(*,'(6x,''socc('',i2,'')='',F18.14)') i, aocc(i)
         enddo
      enddo

!     do i=1,n
!        write(*,*) sum(socc(i,1:nsh))
!     enddo
!    

      end

      subroutine getnsh(j,nsh)
      integer j,nsh
      nsh = 5
      if(j.le.2) nsh=3
      if(j.ge.11.and.j.le.12) nsh=6
      if(j.ge.19.and.j.le.20) nsh=6
      if(j.ge.37.and.j.le.38) nsh=6
      if(j.ge.55.and.j.le.56) nsh=6
      if(j.ge.21.and.j.le.30) nsh=7
      if(j.ge.21.and.j.le.30) nsh=7
      if(j.ge.39.and.j.le.48) nsh=7
      if(j.ge.57.and.j.le.80) nsh=7
      end
