c output the basic molecule info
      subroutine molecho(io,n,nel,na,nb,xyz,rab,at,z)
      use com
      implicit none
      integer n,nel,na,nb
      integer at(n)
      integer io
      real*8  z (n)
      real*8 xyz(3,n)
      real*8 rab(n*(n+1)/2) 

      integer i,j,lin
      character*2 asym
      real*8 cn(n),r,tmp,arg,rcovij

!     standard CN for output
      cn = 0
      do i = 2, n
      do j = 1, i-1 
         r = rab(lin(i,j))
         rcovij=(4./3.)*(rcov(at(i))+rcov(at(j)))
         arg = (r-rcovij)/rcovij
         tmp = 0.5d0 * (1d0 + erf(-7.5d0*arg)) 
         cn(i) = cn(i) + tmp      
         cn(j) = cn(j) + tmp 
      enddo
      enddo

      write(io,*)
      write(io,'(''==============================================='')') 
      write(io,*)'               Molecule Info'
      write(io,'(''==============================================='')') 
      write(io,'('' Number of atoms              : '',i3)') n
      write(io,'('' Number of electrons          : '',i3)') nel
      write(io,'('' Number of alpha electrons    : '',i3)') na
      write(io,'('' Number of beta  electrons    : '',i3)') nb
      write(io,'(''==============================================='')') 
      write(io,*)
      write(io,*)'                                          Coordinates'
      write(io,'(
     .''   Z      CN'',
     .''       x           y          z'')')
      do i=1,n
      write(io,'(a2,i2,f8.3,3f12.6)')
     .asym(at(i)),idint(z(i)),cn(i),xyz(1:3,i)
      enddo   

      end

c output the computed energies             
      subroutine eecho(io,etot,eel,enuc,edisp,eref,exc,ecoul)
      implicit none
      real*8 etot,eel,enuc,edisp,eref,exc,ecoul
      integer io

      integer i
      character*2 asym
      real*8 vir

      write(io,*)
      write(io,'(''==============================================='')') 
      write(io,*)'               Energies'
      write(io,'(''==============================================='')') 
      write(io,'('' repulsion                : '',F18.8)') enuc 
      write(io,'('' dispersion               : '',F18.8)') edisp
      write(io,'('' atomic EXC/core terms    : '',F18.8)') exc  
      write(io,'('' electronic               : '',F18.8)') eel  
      write(io,'('' Coulomb                  : '',F18.8)') ecoul
      if(eref.ne.0)then
      write(io,'('' reference energy         : '',F18.8)') eref 
      endif
      write(io,'('' total                    : '',F18.8)') etot 
      write(io,'(''==============================================='')') 

      end

