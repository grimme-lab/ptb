      character*255 atmp
      real*8 xx(10)
      integer nn

      open(unit=1,file='.QREF')
  10  read(1,'(a)',end=100) atmp
      call readl(atmp,xx,nn)
      if(nn.eq.6.and.xx(1).gt.1.d-3) call system('touch alpha')
      goto 10

100   continue
      close(1)
      end
