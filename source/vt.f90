!--------------------------------------------------------
! compute atomic/shell energy increments/VTC corrections
!--------------------------------------------------------

subroutine gtbenergy(pr,n,nao,at,z,xyz,rab,pa,psh,S,P,eref,etot)
      use bascom
      use parcom
      use com
      use dftd4     
      implicit none
      integer n,nao
      logical,intent(in)     :: pr          
      integer,intent(in)     :: at(n)       
      real*8, intent(in)     :: z(n)
      real*8, intent(in)     :: xyz(3,n)
      real*8, intent(in)     :: rab(n*(n+1)/2) 
      real*8, intent(in)     :: pa(n)                     ! atom populations
      real*8, intent(in)     :: psh(10,n)                 ! shell populations
      real*8, intent(in)     :: S(nao*(nao+1)/2)          ! overlap integrals 
      real*8, intent(in)     :: P(nao*(nao+1)/2)          ! P 
      real*8, intent(in)     :: eref                      ! ref. energy  
      real*8, intent(out)    :: etot                      ! total energy 

      real*8, parameter      :: erfs = -2.0d0             ! CN
      integer i,j,k,ij,ati,atj,ia,ib,ish,jsh
      real*8 cn(n),scal(10,n)
      real*8 glob8,glob9,glob10,glob11,glob12
      real*8 zi, zj, ai, aj, zij
      real*8 r0i,r0j,aa,aa2,rr,r0ab
      real*8 qi,qj,hi,hj,hij
      real*8 s8,a1,a2
      real*8 evt,eel,edisp
      real*8,allocatable :: T(:),SS(:),xnorm(:)

      allocate(SS(nao*(nao+1)/2),T(nao*(nao+1)/2),xnorm(nao))     

      call ncoord_erf(n,at,rab,erfs,cn)

      call sscal(n,at,psh,glob_par(19),glob_par(20),scal)
      call modbasd(n,at,scal)                             ! scale exponents shell/atom-wise 
      call stint(n,nao,at,xyz,rab,SS,T,xnorm)             ! scaled T, SS dummy
      call modbas(n,at,4)                                 ! back


!     force correction to VT (nuclear repulsion analogon)
      glob9 =glob_par(15)
      glob10=glob_par(16)
      glob11=glob_par(17)
      glob12=glob_par(18)
      k=0
      evt=0
      do i=1,n
         ati=at(i)
         qi =z(i)-pa(i)
         zi = ener_par1(1,ati) 
         ai =(ener_par1(2,ati)+ener_par1(4,ati)*cn(i))*(1d0+glob10*qi)
         r0i=(ener_par1(3,ati)+ener_par1(5,ati)*cn(i))*(1d0+glob9 *qi)
         do j=1,i-1
            k  = k + 1
            atj=at(j)
            qj =z(j)-pa(j)
            zj = ener_par1(1,atj) 
            aj =(ener_par1(2,atj)+ener_par1(4,atj)*cn(j))*(1d0+glob10*qj)
            r0j=(ener_par1(3,atj)+ener_par1(5,atj)*cn(j))*(1d0+glob9 *qj)
            aa    = sqrt(ai * aj)
            r0ab  = (r0i + r0j)  !*(1d0+qi*qj*glob8) !+qi*qi*qj*qj*glob11)
            rr    = r0ab / rab(k)
            aa2   = aa*0.5d0         
            zij   = sqrt(zi*zj)*(1d0+qi*qj*glob11+qi*qi*qj*qj*glob12)
            evt   = evt  + zij*(rr**aa-rr**aa2) ! general LJ
         enddo
         k = k +1
      enddo
 
      ij = 0
      do i=1,nao 
         ati = at(aoat(i))
         ish = shell2ao(i)
         hi  = ener_par2(ish,ati)
         do j=1,i  
            ij  = ij + 1       
            atj = at(aoat(j))
            jsh = shell2ao(j)
            hj  = ener_par2(jsh,atj)
            hij = 0.5d0*(hi + hj)
            T(ij) = T(ij) + hij * S(ij)                                   
         enddo
      enddo

      call energy(nao,T,P,eel)  ! E_T

!     dispersion
      s8 = 0.5d0   !glob_par(8)
      a1 = 0.45d0  !glob_par(9)
      a2 = 4.6d0   !glob_par(10)
      cn = z - pa  ! q
      call dftd4_dispersion(at, xyz, cn, 1.0d0, s8, 1.0d0, a1, a2, 3.0d0, 2.0d0, edisp)
   
!     add up scaled T + VT force + disp
      etot = -eel*(1.3d0+glob_par(14)) + evt + edisp  

      if(pr) then
         write(*,*) 'dispersion              :',edisp
         write(*,*) 'atomic VT force         :',evt 
         write(*,*) 'kinetic energy          :',eel
         write(*,*) 'total energy (calc/ref) :',etot,eref
         if(abs(eref).gt.1.d-9) &
         write(*,*) 'energy error per /atom  :',(etot-eref)/float(n)
      endif

      end

!! ------------------------------------------------------------------------

