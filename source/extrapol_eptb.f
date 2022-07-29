      implicit none   
      real*8 e0, e1, e
      real*8 fia0  
      real*8 p(1:10)
      real*8 c,fe,ff,ex
    
      open(unit=1,file='.tmpx')
      read(1,*) e0
      read(1,*) e1
      read(1,*) fia0
      close(1)

      open(unit=2,file='~/.atompara')
      read(2,*) p(1:10)
      close(2)

      c=100.0*(e1-e0)/(0.5*(e1+e0)) ! % E change

      fe = 1.0/(1.0+(   c/p(5))**16)  ! % E change damping
      ff = 1.0/(1.0+(fia0/p(6))**16)  ! max fia damping

!        E part   fia part  base
      ex=p(7)*fe + p(8)*ff + p(1)  ! smaller p(1) (i.e. 0) is better for MG, larger for TMs

      e = e1 + (e1-e0) * ex

!     write(*,*) 'fe ff',fe,ff
      write(*,'(4F20.8)') e, c, fia0, ex
      end 

