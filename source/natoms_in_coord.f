      real*8  ,allocatable :: xyz(:,:)
      integer, allocatable :: at(:)
      integer n
      character*80 fname

      call getarg(1,fname)
      call rd0(fname,n)
      allocate(at(n),xyz(3,n))
      call rd(.false.,fname,n,xyz,at)

      write(*,*) n 

      end
