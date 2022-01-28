!--------------------------------------------------------
! compute atomic/shell energy increments/VTC corrections
!--------------------------------------------------------

subroutine gtbenergy(pr,n,nao,at,z,xyz,rab,pa,psh,P,wbo,eref,etot)
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
      real*8, intent(in)     :: P(nao*(nao+1)/2)          ! P 
      real*8, intent(in)     :: wbo(n,n)                  ! WBO 
      real*8, intent(in)     :: eref                      ! ref. energy  
      real*8, intent(out)    :: etot                      ! total energy 

      real*8, parameter      :: erfs    = -2.0d0             ! CN
!     real*8, parameter      :: gamscal = -1.00d0            ! 
      integer i,j,k,ij,ati,atj,ia,ib,ish,jsh,ksh,ii,jj,lin
      real*8 cn(n),scal(10,n)
      real*8 zi, zj, ai, aj, zij, ff, pol, gamscal
      real*8 r0i,r0j,r2, r0ab,damp
      real*8 qi,qj,hi,hj,hij,gish,gjsh,xk,tmp
      real*8 s8,a1,a2
      real*8 atocc(10),qsh(nao)
      real*8 enuc,edisp,e1,ecoul,eel
      real*8,allocatable :: gab(:,:)
      real*8,allocatable :: SS(:),xnorm(:)

      call ncoord_erf(n,at,rab,erfs,cn)

      gamscal=glob_par(20)
!!!!!!!!!!!!!!
! dispersion
!!!!!!!!!!!!!!
      s8 = 0.5d0   !glob_par(8)
      a1 = 0.45d0  !glob_par(9)
      a2 = 4.6d0   !glob_par(10)
      cn = z - pa  ! q
      call dftd4_dispersion(at, xyz, cn, 1.0d0, s8, 1.0d0, a1, a2, 3.0d0, 2.0d0, edisp)

      e1 = 0 
      ksh= 0
      do i=1,n
         ati= at(i)
         qi = z(i)-pa(i)
         ff = (pa(i) / z(i))**0.5d0        ! atom specific energy increment normalization factor (makes E(H+)=0)
         e1 = e1 + ff* ener_par2(1,ati)    ! free atom
         e1 = e1 + ff* ener_par2(2,ati)*qi ! charge correction
         call shellocc_ref(ati,atocc)      ! ref. atomic pop.
         do ish=1,bas_nsh(ati)
            ksh = ksh + 1
            qsh(ksh) = atocc(ish)-psh(ish,i)
         enddo
      enddo

!!!!!!!!!!!!!!!!!!
! Coulomb energy
!!!!!!!!!!!!!!!!!!
      allocate(gab(nsh,nsh))
      ii = 0
      do i=1, n
         ati = at(i)
         qi  = z(i)-pa(i)
         do ish=1, bas_nsh(ati)
         ii = ii + 1
         gish = gam(ati) * shell_cnf3(ish,ati) * (1d0 + gamscal*qi) 
         jj = 0
         do j=1,n
            k = lin(j,i)
            r2= rab(k)**2
            atj = at(j)
            qj  = z(j)-pa(j)
            do jsh=1, bas_nsh(atj)
               jj = jj + 1
               if (jj.gt.ii) cycle
               gjsh = gam(atj) * shell_cnf3(jsh,atj) * (1d0 + gamscal*qj) 
               xk= 2d0 /(1d0/gish + 1d0/gjsh)    ! harm. av.
               gab(jj,ii)=1d0/sqrt(r2+1d0/xk**2) ! Ohno-Klopman-Mataga average
               gab(ii,jj)=gab(jj,ii)
            enddo
         enddo
         enddo
      enddo
      ecoul=0
      do i=1,nsh-1
         do j=i+1,nsh
            ecoul =ecoul + qsh(i)*qsh(j)*gab(j,i)
         enddo
      enddo
      ecoul=ecoul*2.0d0
      do i=1,nsh
         ecoul =ecoul + qsh(i)*qsh(i)*gab(i,i)
      enddo
      ecoul=ecoul*0.5d0

!!!!!!!!!!!!!!
! electronic
!!!!!!!!!!!!!!

      allocate(SS(nao*(nao+1)/2),xnorm(nao))

      call sscal(n,at,psh,glob_par(15),glob_par(16),scal) ! 15=0.0336328354 16=-0.0085599299 17=0.6572780940 18=-0.0017648814 19=0.8565958171 20=-1.3709605352
      call modbasd(n,at,scal)                             ! scale exponents shell/atom-wise 
      call sint(n,nao,at,xyz,rab,SS,xnorm)                ! scaled T, SS dummy
      call modbas(n,at,4)                                 ! back

      ij = 0
      eel= 0
      do i=1,nao  
         ia = aoat(i)
         ati= at(ia)
         ish= shell2ao(i)
         hi = ener_par1(ish,ati)
         do j=1,i-1
            ij = ij + 1
            ib = aoat(j)
            atj= at(ib)
            jsh= shell2ao(j)
            hj = ener_par1(jsh,atj)
            hij= hi+hj
            pol = ((hi-hj)/hij)**2
            if(ia.ne.ib) then            ! different atoms
               tmp = hij * (1d0-pol*glob_par(18))
            else                         ! same atoms
               tmp = hij 
               if(ish.ne.jsh) then       ! s-s', p-p' off-diagonal (there is no d-d' in MG VDZP)
                  tmp = hij * glob_par(19)
               endif
            endif
!           H(ij) = tmp * SS(ij) 
            eel = eel + tmp * SS(ij) * P(ij)
         enddo
         ij = ij + 1
!        H(ij) = 2d0 * hi * SS(ij) 
         eel = eel + hi * SS(ij) * P(ij)
      enddo
      eel = eel * 2d0
!     call energy(nao,H,P,eel) 

!!!!!!!!!!!!!!
! repulsion
!!!!!!!!!!!!!!
      enuc=0
      k   =0
      do i=1,n
         ati=at(i)
         qi = z(i)-pa(i)
         zi = ener_par2(4,ati) * (1d0+ener_par2(5,ati)*qi)  ! q has good effect
         r0i= ener_par2(3,ati) 
         do j=1,i-1
            k   = k + 1 
            atj = at(j)
            qj  = z(j)-pa(j)
            zj  = ener_par2(4,atj) * (1d0+ener_par2(5,atj)*qj)
            r0j = ener_par2(3,atj) 
            r0ab= (r0i + r0j) *(1d0+qi*qj*glob_par(14)) !=0.0907067762 
            damp= erf(-r0ab*rab(k))+1d0
            zij = zi * zj * (1d0+glob_par(17)*wbo(j,i))
            enuc= enuc + damp*zij/rab(k)
         enddo
         k = k +1
      enddo


!     add up 
      etot = eel + edisp + e1 + ecoul + enuc

      if(pr) then
         write(*,*) 'dispersion              :',edisp
         write(*,*) 'atomic increment        :',e1  
         write(*,*) 'nuclear repulsion       :',enuc 
         write(*,*) 'Coulomb                 :',ecoul
         write(*,*) 'electronic energy       :',eel
         write(*,*) 'total energy (calc/ref) :',etot,eref
         if(abs(eref).gt.1.d-9) &
         write(*,*) 'energy error per /atom  :',(etot-eref)/float(n)
      endif

      end

