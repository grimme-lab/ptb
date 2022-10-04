!--------------------------------------------------------
! compute total energy bei TB2 type expression and PTB
! input, SG 6/22
!        read(1,*) ener_par1 (1:10,j)    ! 1 -10       
!        read(1,*) ener_par2 (1:10,j)    ! 11-20       
!        read(1,*) expscal  (1,1:10,j)   ! 21-30  
!        read(1,*) ener_par6 (1:10,j)    ! 31-40  
!        read(1,*) ener_par4 (1:10,j)    ! 41-50       
!        read(1,*) ener_par5 (1:10,j)    ! 51-60       
!--------------------------------------------------------

subroutine ptb_energy(pr,n,nao,nopen,at,z,xyz,rab,pa,psh,wbo,S,P,eref,etot)
      use parcom
      use com
      use bascom
      use dftd4 
      use aescom, only: qm,dipm,qp
      implicit none
      integer n,nao,nopen
      logical,intent(in)     :: pr          
      integer,intent(in)     :: at(n)       
      real*8, intent(in)     :: z(n)
      real*8, intent(in)     :: xyz(3,n)
      real*8, intent(in)     :: rab(n*(n+1)/2) 
      real*8, intent(inout)  :: pa(n)                     ! atom populations
      real*8, intent(inout)  :: psh(10,n)                 ! shell populations
      real*8, intent(inout)  :: wbo(n,n)                  ! WBO (not used)
      real*8, intent(in)     :: S(nao*(nao+1)/2)          ! exact S     
      real*8, intent(inout)  :: P(nao*(nao+1)/2)          ! PTB density
      real*8, intent(out)    :: eref,etot                 ! total energy 

      real*8,parameter   :: eTemp  =300.0d0   ! electronic temp. 
      integer i,j,k,l,ij,kl,ati,atj,ia,ib,ish,jsh,ii,jj,lin,li,lj
      integer llao(4)
      data llao / 1,3,5,7 /
      integer, parameter :: idx(3, 3) = reshape([1, 2, 4, 2, 3, 5, 4, 5, 6], [3, 3])
      real*8 ea,erep,eel,edisp,ecoul,exc,ewbo,eq3,eaes,eac,ecp,ge
      real*8 r0i,r0j,r0ab,damp
      real*8 zi,zj,ff,xk,r,t8,t9,arg,rcovij,ten(3)
      real*8 hi,hj,hij,pol,tmp,tmp2,tmp3,keav,ssh,nel
      real*8 s8,a1,a2,s9
      real*8 t0,t1,w0,w1
      real*8 xx(10),par(5)
      character*255 atmp
      logical ex
      real*8,allocatable :: SS(:),Hdiag(:),xnorm(:),scal(:,:),spsh(:,:),q(:)
      real*8,allocatable :: sbo(:),cn(:),cn1(:),cn2(:),eps(:)


      allocate(Hdiag(nao),cn1(n),cn(n),cn2(n),sbo(n),q(n),eps(nao), &
     &         SS(nao*(nao+1)/2),xnorm(nao),scal(10,nsh),spsh(10,n), &
     &         qm(n),dipm(3,n),qp(6,n))

      eref =0
      edisp=0
      ex   =.false.

      if(pr)then
      inquire(file='control', exist=ex)
      if(ex) then  
      open(unit=10,file='control')
25    read(10,'(a)',end=30) atmp
      if(index(atmp,'subenergy').ne.0) then
      read(10,'(a)',end=30) atmp
      call readl(atmp,xx,ii)
      if(ii.ge.3) eref=xx(1)-xx(7)  ! remove Edisp
      goto 30
      endif
      goto 25
30    close(10)
      endif
      endif

      spsh = 0
      open(unit=11,file='ptb_dump',form='unformatted')
      read(11) pa  
      read(11) psh
      read(11) wbo
      read(11) P  
      read(11) qm 
      read(11) dipm
      read(11) qp
      read(11) eps  
      read(11) xnorm
      if(nopen.ne.0) read(11) spsh
      close(11)
 
!!!!!!!!!!!!!!
! dispersion
!!!!!!!!!!!!!!

      q = z - pa ! q from PTB pop  

!                                      s8  s9           a1     a2    beta_1/2 (orig 3.0,2.0)
!     call dftd4_dispersion(at,xyz,q,1.0d0,0d0,1d0,0.2464d0,4.4737d0,3d0,2d0,edisp)  ! wB97X-3c values

      s8 = 1.00d0  !glob_par(7)
      a1 = 0.55d0  !glob_par(8)
      a2 = 4.70d0  !glob_par(9)
      s9 = 2.0d0
      call dftd4_dispersion(at,xyz,q,1d0,s8,s9,a1,a2,3d0,2d0,edisp)  ! fitted values

      if(ex) eref = eref + edisp ! take same disp for fit (i.e. not considered in fit but in energy output)
                                 ! the fit is thus to wB97X/vDZP for energies

!!!!!!!!!!!!!!
! CNs         
!!!!!!!!!!!!!!

      call ncoord_erf(n,at,rab,-2d0,cn1) ! all erf damping argument are -2!

      cn = 0
      do i = 2, n
      do j = 1, i-1 
         r = rab(lin(i,j))
         rcovij=abs(shell_cnf4(3,at(i)))+abs(shell_cnf4(3,at(j))) ! PTB values
         arg = (r-rcovij)/rcovij
         tmp = 0.5d0 * (1d0 + erf(-2d0*arg)) 
         cn(i) = cn(i) + tmp      
         cn(j) = cn(j) + tmp 
      enddo
      enddo

!!!!!!!!!!!!!!
! atomic   
!!!!!!!!!!!!!!

      do i=1,n
         wbo(i,i)=0
         tmp2=0
         do j=1,n
            tmp2=tmp2+wbo(j,i)
         enddo
         sbo(i)=tmp2
      enddo

      ea = 0 
      eq3= 0 
      eac= 0
      do i=1,n
         ati = at(i)
         if(abs(ener_par1(2,ati)).lt.1.d-6) stop 'parameter missing'

         ea = ea +   ener_par1(1,ati)*z(i)*(1d0+ener_par1(4,ati)*cn1(i)) ! DFT energy shift

         eq3 = eq3 + ener_par1(2,ati)*q(i)**3   &               ! third order diag as in GFN2
     &             + ener_par1(3,ati)*q(i)**4   &               ! new
     &             + glob_par(18)    *q(i)**5                   ! new
         t8 = 0
         do k = 1,3
         do l = 1,3
            kl = idx(l,k)
            t8 = t8 + qp(kl,i)*qp(kl,i)
         enddo
         enddo
         t9 =sum(dipm(1:3,i)**2)
         eac = eac + ener_par5(9,ati)*t9 + ener_par5(10,ati)*t8 ! AXC as in GFN2

!        tmp=0
!        do ish=1,bas_nsh(ati)
!           tmp = tmp + spsh(ish,i)
!        enddo
!        esp = esp + ener_par5(7,ati)*tmp    ! spin density correction

      enddo

!!!!!!!!!!!!!!
! Coulomb        
!!!!!!!!!!!!!!

      call shell_es(n,at,rab,q,psh,ecoul)

      call aniso_es(n,at,rab,xyz,eaes) 
       
!!!!!!!!!!!!!!
! ECP
!!!!!!!!!!!!!!
      
      call eecp(n,nao,at,xyz,rab,xnorm,P,ecp)

!!!!!!!!!!!!!!
! electronic
!!!!!!!!!!!!!!

!                                         call timing(t0,w0)           
      call shscalE(n,at,psh,scal)
      call modbasd(n,at,scal)              ! scale exponents shell-wise with psh dep.
      call sint(n,nao,at,xyz,rab,SS,xnorm) 
      call modbas(n,at,4) 
!                                         call timing(t1,w1)           
!                                         call prtime(6,t1-t0,w1-w0,'S')

!     atomic H0 
      nel = 0
      ii  = 0
      do i = 1, n
         ati = at(i)
         do ish=1,bas_nsh(ati)
            tmp = ener_par2(ish,ati) + ener_par4(ish,ati)*(cn1(i)+cn(i)*ener_par1(9,ati)) ! similar (but not identical) to PTB 
            nel = nel + psh(ish,i)
            do j=1,llao(bas_lsh(ish,ati)+1) ! AO loop
               ii = ii + 1
               Hdiag(ii) = tmp
            enddo
         enddo
      enddo


!     H0
      eel= 0
      ij = 0
      do i=1,nao 
      ia = aoat(i)
      ati= at(ia)
      ish= shell2ao(i)
      li = bas_lsh(ish,ati)
      hi = Hdiag(i)
      do j=1,i  
         ij = ij + 1
         ib = aoat(j)
         hj = Hdiag(j)
         hij= hi+hj
         ssh= hij * SS(ij)
         atj= at(ib)
         jsh= shell2ao(j)
         lj = bas_lsh(jsh,atj)
         if(ia.ne.ib) then            ! different atoms
            pol = (hi-hj)**2        ! denominator removed (unstable in fit)
            keav= sqrt(ener_par5(li+1,ati) * ener_par5(lj+1,atj)) 
            tmp = ssh * keav * (1d0-pol*glob_par(19))   ! optimized, R dep. removed 945 ! tmp = tmp * (1d0+0.001d0/rab(lin(ia,ib)))
         else                         ! same atoms
            if(ish.ne.jsh) then       ! s-s', p-p', d-d' off-diagonal, li=lj because S=0 otherwise
               tmp = ssh * ener_par2(8+li,ati) * (1d0+sqrt(cn1(ia))*ener_par1(6,ati))  ! new CN dep.
            else
               tmp = ssh 
            endif
         endif
         if(i.ne.j) then
            eel = eel + tmp * P(ij) * 2d0
         else
            eel = eel + tmp * P(ij)
         endif
      enddo
      enddo

!!!!!!!!!!!!!!
!XC
!!!!!!!!!!!!!!

!                                             call timing(t0,w0)           
   call shscalE2(n,at,scal)
   call modbasd(n,at,scal)       
   call sint(n,nao,at,xyz,rab,SS,xnorm) 
   call modbas(n,at,4) 
   call epauli2(n,nao,at,psh,SS,P,cn1,q,Hdiag,exc) ! add valence X correction based on three-index ECP formula 
                                                   ! three parameters per element (one in shscalE2)
!                                             call timing(t1,w1)           
!                                             call prtime(6,t1-t0,w1-w0,'X')


!!!!!!!!!!!!!!
! repulsion
! and WBO term
!!!!!!!!!!!!!!
      erep=0
      ewbo=0
      k   =0
      do i=1,n
         ati=at(i)
         r0i= ener_par1(5,ati) + cn1(i)*ener_par1(8,ati)*2.0
         zi = z(i) * (1d0+ener_par1(7,ati)*(q(i)-ener_par5(4,ati)*q(i)**2))  ! important charge correction
!        hi = 0.05d0
!        if(ati.le.2) hi = 1.0d0
         do j=1,i-1
            k   = k + 1 
            atj = at(j)
            hj  = 0.05d0
!           if(atj.le.2) hj = 1.0d0
!           ssh = 0.5d0*(hi+hj)
            ssh = 0.5d0*(expscal(1,8,ati) + expscal(1,8,atj))

            r0j = ener_par1(5,atj) + cn1(j)*ener_par1(8,atj)*2.0
            zj  = z(j) * (1d0+ener_par1(7,atj)*(q(j)-ener_par5(4,atj)*q(j)**2))

            r0ab= sqrt(r0i * r0j)                                      ! slightly better than amean
            damp= 0.5d0*(erf(-2d0*(rab(k)-r0ab)/r0ab**2)+1d0) 
            t8  = abs(wbo(j,i))**ssh            
            tmp = 1d0 +(ener_par4(8,atj)+ener_par4(8,ati))*t8

            erep= erep + damp * zi * zj * tmp / rab(k)

            tmp = 0.5d0*(ener_par5(7,ati)+ener_par5(7,atj))
            tmp2= sqrt(ener_par5(5,ati)*ener_par5(5,atj))
            ewbo= ewbo + tmp2                                &    ! this is a less empirical version of the "+U" in PTB
     &                 * (wbo(j,i)+tmp*wbo(j,i)**2)               ! here its energy and WBO is (PS)^2 

!           write(*,*) (ener_par5(5,atj)+ener_par5(5,ati))*wbo(j,i),(ener_par5(5,atj)+ener_par5(5,ati))*tmp*wbo(j,i)**2
                  
         enddo
         k = k + 1
      enddo

!!!!!!!!!!!!!!
! add up 
!!!!!!!!!!!!!!

!     ii=idint(nel)
!     call free_en(nao,ii,nopen,etemp,eps,ge)

      eel = eel * 0.4d0           
      erep=erep * glob_par(20)
      etot = eel + erep + ea + ecoul + eaes + eq3 + eac + exc + ewbo + edisp + ecp 

      if(pr) then    ! some output of energies and atomic descriptors
         write(*,*) 
         write(*,*) '             -----------------------'
         write(*,*) '             |  EPTB energy model  |'
         write(*,*) '             |       SG 6/22       |'
         write(*,*) '             -----------------------'
         write(*,'(''    atom             q            CN        sum BO  dipole r^2 aniso  spin pop'')')
         do i=1,n
            tmp=0
            do ish=1,bas_nsh(at(i))
               tmp = tmp + spsh(ish,i)
            enddo
            t8 = 0
            do k = 1,3
            do l = 1,3
            kl = idx(l,k)
            t8 = t8 + qp(kl,i)*qp(kl,i)
            enddo
            enddo
            t9=sum(dipm(1:3,i)**2)
            write(*,'(2i5,f5.1,3x,11f8.3)') i,at(i),z(i),z(i)-pa(i),cn1(i),cn(i),sbo(i),t9,t8,tmp             
         enddo
         write(*,*) 
         write(*,'(''# electrons             :'',2F16.8)')nel   
!        write(*,'(''free energy at T        :'',2F16.8)')ge,etemp
         write(*,*) 
         write(*,'(''  electronic (H0*P)     :'',2F16.8)')eel   
         write(*,'(''  WBO                   :'',2F16.8)')ewbo  
         write(*,'(''electronic+WBO          :'',2F16.8)')ewbo+eel
         write(*,'(''  Coulomb 2nd-order     :'',2F16.8)')ecoul
         write(*,'(''  anisotropic ES        :'',2F16.8)')eaes 
         write(*,'(''  third-order ES        :'',2F16.8)')eq3  
         write(*,'(''total ES                :'',2F16.8)')eq3 + eaes + ecoul
         write(*,'(''  XC (non-local)        :'',2F16.8)')exc  
         write(*,'(''  multipole XC          :'',2F16.8)')eac 
         write(*,'(''  atomic                :'',2F16.8)')ea  
         write(*,'(''total XC + atomic       :'',2F16.8)')ea+exc+eac
         write(*,'(''ECP                     :'',2F16.8)')ecp  
         write(*,'(''dispersion              :'',2F16.8)')edisp
!        write(*,'(''spin density            :'',2F16.8)')esp
         write(*,'(''nuclear repulsion       :'',2F16.8)')erep 
         write(*,*) 
         if(abs(eref).gt.1.d-9) then
         write(*,'(''total energy (calc/ref) :'',2F16.8)')etot,eref
         write(*,'(''energy error per /atom  :'',2F16.8)')(etot-eref)/float(n)
         else
         write(*,'(''total energy            :'',2F16.8)')etot
         endif
      endif

      deallocate(qm,dipm,qp)

      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! second order Coulmob energy with special gammas
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine shell_es(n,at,rab,q,qsh,ecoul)
   use iso_fortran_env, only : wp => real64
   use bascom
      implicit none 
      integer, intent(in)  :: n
      integer, intent(in)  :: at(n)
      real*8,  intent(in)  :: rab(n*(n+1)/2)  
      real*8,  intent(in)  :: q(n)
      real*8,  intent(in)  :: qsh(10,n)
      real*8,  intent(out) :: ecoul    

      integer i,j,ati,ish,iish
      real*8  atocc(10)
      real*8  qshtmp(nsh),gab(nsh,nsh)

      call setgab3(n,at,rab,q,gab) ! PTB except for other par(3) (=0 in PTB)

      iish = 0
      do i=1,n
         ati = at(i)
         do ish=1,bas_nsh(ati)
            iish = iish + 1
            call shellocc_ref(ati,atocc) ! ref. atomic pop.
            qshtmp(iish)=atocc(ish)-qsh(ish,i)
         enddo
      enddo

      ecoul=0
      do i=1,iish-1
         do j=i+1,iish
            ecoul =ecoul + qshtmp(i)*qshtmp(j)*gab(j,i)
         enddo
      enddo
      ecoul=ecoul*2.0d0
      do i=1,iish
         ecoul =ecoul + qshtmp(i)*qshtmp(i)*gab(i,i)
      enddo
      ecoul=ecoul*0.5d0

end  

!! ------------------------------------------------------------------------
!! ------------------------------------------------------------------------

subroutine eecp(n,nao,at,xyz,rab,norm,P,e)
      use cbascom
      use  bascom
      use  parcom
      use gtb_la, only : la_gemm
      implicit none          
      integer, intent(in)   :: nao,n,at(n)
      real*8,  intent(in)   :: xyz(3,n)
      real*8,  intent(in)   :: rab(n*(n+1)/2) 
      real*8,  intent(in)   :: norm(nao)    
      real*8,  intent(in)   :: P(nao*(nao+1)/2)    
      real*8,  intent(out)  :: e

      integer i,j,k,l,m,nl,nn,atn,jsh,llao2(0:3),ia,ib
      data llao2/1,3,5,7 /
      real*8 vecp, ddot
      real*8,allocatable :: Scv(:,:), stmp(:,:), xtmp(:,:)

      e = 0

      if(cnsao.eq.0) return

      allocate(Scv(cnsao,nao),stmp(cnsao,nao),xtmp(nao,nao))

      call csint(n,nao,at,xyz,rab,norm,Scv) ! core val overlap ints

!     N^2 step
      do i=1,nao                          
         m=0 
         do nl=1,ncorelist                     ! all atoms with core
            nn=corelist(nl)
            atn=at(nn)
            do jsh=1,cbas_nsh(atn)             ! core shells of atom nn
               do l=1,llao2(cbas_lsh(jsh,atn)) ! AOs of core shell jsh
                  m = m + 1
                  stmp(m,i)= clev(jsh,atn) * Scv(m,i) * expscal(1,10,atn) ! * shell_cnf1(10,atn)
               enddo
            enddo
         enddo
      enddo

!     N^3 step
      call la_gemm('T','N',nao,nao,cnsao,1.0d0,Scv,cnsao,stmp,cnsao,0.0d0,xtmp,nao)

      k = 0 
      do i=1, nao 
         do j=1, i-1
            k = k + 1 
            e = e + xtmp(j,i) * P(k)
         enddo
         k = k + 1
         e = e + xtmp(i,i) * P(k) * 0.5d0
      enddo

      e = 2d0 * e

end

!! ------------------------------------------------------------------------
!  add Pauli term to H from 2. iter
!! ------------------------------------------------------------------------

subroutine epauli2(n,nao,at,psh,S,P,cn,q,Hdiag,e)
      use  bascom
      use  parcom
      use gtb_la, only : la_symm
      implicit none          
      integer, intent(in)   :: nao,n,at(n)
      real*8,  intent(in)   :: psh(10,n)
      real*8,  intent(in)   :: S(nao*(nao+1)/2)    
      real*8,  intent(in)   :: P(nao*(nao+1)/2)    
      real*8,  intent(in)   :: Hdiag(nao)    
      real*8,  intent(in)   :: cn(n)    
      real*8,  intent(in)   :: q(n)    
      real*8,  intent(out)  :: e

      integer i,j,k,l,m,nl,atk,jsh,llao2(0:3)
      data llao2/ 1,3,5,7 /
      real*8 f1,dum
      real*4,allocatable :: stmp(:,:), sdum(:,:), xtmp(:,:)

      allocate(stmp(nao,nao),sdum(nao,nao),xtmp(nao,nao))

      call blowsym84(nao,S,sdum)

!     N^2 step
      do i=1,nao                          
         m=0 
         do k=1,n       
            atk=at(k)
            dum=ener_par1(10,atk)+cn(k)*ener_par5(6,atk)
            do jsh=1,bas_nsh(atk)          
               l =bas_lsh(jsh,atk)
               nl=llao2(l)
               f1=psh(jsh,k) * dum / dble(nl)
               do l=1,nl                  
                  m = m + 1
                  stmp(m,i)= Hdiag(m) * sdum(m,i) * f1
               enddo
            enddo
         enddo
      enddo

!     N^3 step
      call la_symm('L','L',nao,nao,1.0e0,sdum,nao,stmp,nao,0.0e0,xtmp,nao)   

      k = 0 
      do i=1, nao 
         do j=1, i-1
            k = k + 1 
            e = e + xtmp(j,i) * P(k)
         enddo
         k = k + 1
         e = e + xtmp(i,i) * P(k) * 0.5d0
      enddo

      e = 2d0 * e

end

!! ------------------------------------------------------------------------
!  scale factors for S'
!! ------------------------------------------------------------------------

subroutine shscalE(n,at,psh,scal)
   use iso_fortran_env, only : wp => real64
   use parcom
   use bascom
   implicit none 
!! ------------------------------------------------------------------------
!  Input
!! ------------------------------------------------------------------------
   integer, intent(in)    :: n                  ! number of atoms 
   integer, intent(in)    :: at(n)              ! ordinal number of atoms
   real(wp),intent(in)    ::psh(10,n)           ! shell populations       
!! ------------------------------------------------------------------------
   real(wp),intent(out)   ::scal(10,n)          ! exponent scaling factors
!! ------------------------------------------------------------------------
   integer  :: i,ish,ati
   real(wp) :: atocc(10), qa

   do i=1,n
      ati = at(i)
      call shellocc_ref(ati,atocc) ! ref. atomic pop.
      do ish=1,bas_nsh(ati)
         qa = atocc(ish)-psh(ish,i)
         scal(ish,i) = expscal(1,ish,ati) * (1d0 + ener_par4(9,ati)*(qa-0.1d0*qa**3)) ! ^3 term ensures that factor stays positive -0.06^3
         if(scal(ish,i).lt.0.05) scal(ish,i)=0.05
         if(scal(ish,i).gt.20.0) scal(ish,i)=20.0
      enddo
   enddo

end

subroutine shscalE2(n,at,scal)
   use iso_fortran_env, only : wp => real64
   use parcom
   use bascom
   implicit none 
!! ------------------------------------------------------------------------
!  Input
!! ------------------------------------------------------------------------
   integer, intent(in)    :: n                  ! number of atoms 
   integer, intent(in)    :: at(n)              ! ordinal number of atoms
!! ------------------------------------------------------------------------
   real(wp),intent(out)   ::scal(10,n)          ! exponent scaling factors
!! ------------------------------------------------------------------------
   integer  :: i,ish,ati,l

   do i=1,n
      ati = at(i)
      do ish=1,bas_nsh(ati)
         l = bas_lsh(ish,ati)
         scal(ish,i) = expscal(1,9,ati) * (1d0-dble(l)*0.00000d0)
      enddo
   enddo

end
!! ------------------------------------------------------------------------
!     DFTB second order term J matrix (other mean than in PTB)
!! ------------------------------------------------------------------------

subroutine setgab3(n,at,rab,q,gab)
   use bascom
   use parcom
   use com
      implicit none 
      integer, intent(in)  :: n
      integer, intent(in)  :: at(n)
      real*8,  intent(in)  :: rab(n*(n+1)/2)  
      real*8,  intent(in)  :: q(n)
      real*8,  intent(out) :: gab(nsh,nsh)

      integer i,j,k,ati,atj,ish,jsh,ii,jj,lin
      real*8 gish,gjsh,xk,r2,geff(n),cok,cmn

!     cok=1.05d0        ! global parameter
!     cmn=1.0d0-cok

!     do i=1,n
!        geff(i) = (1d0 - gsc*(q(i)-0.00*q(i)**2))*gam(at(i))
!     enddo

      ii = 0
      do i=1, n
      ati = at(i)
      do ish=1, bas_nsh(ati)
         ii = ii + 1
         gish = ener_par6(ish,ati) * gam(ati)
         jj = 0
         do j=1,n
            k = lin(j,i)
            r2= rab(k)**2
            atj = at(j)
            do jsh=1, bas_nsh(atj)
               jj = jj + 1
               if (jj.gt.ii) cycle
               gjsh = ener_par6(jsh,atj) * gam(atj)
!              if(i.eq.j)then
               xk   = sqrt  (gish * gjsh) 
!              else
!              xk   = 0.5d0*(gish + gjsh) 
!              xk   = 2d0/(1d0/gish + 1d0/gjsh) 
!              endif
               gab(jj,ii)= 1d0/sqrt(r2+1d0/(xk**2 + 1d-6))
!              gab(jj,ii)= cok/sqrt(r2+1d0/(xk**2+1d-6)) + cmn/(rab(k)+1d0/(xk+1d-6))
               gab(ii,jj)= gab(jj,ii)
            enddo
         enddo
      enddo
      enddo

end

subroutine tensav(a,e)
      use gtb_la, only : la_syev
      real*8 a(6),e(3)
      integer info,lwork
      real*8, allocatable ::aux(:),vecs(:,:)

      lwork=37
      allocate (vecs(3,3),aux(lwork))

      k=0
      do i=1,3
         do j=1,i
            k=k+1
            vecs(j,i)=a(k)
            vecs(i,j)=a(k)
         enddo
      enddo
      call la_syev ('V','U',3,vecs,3,e,aux,lwork,info)
end

subroutine free_en(ndim,nel,nopen,et,eps,g)
      implicit none
      integer ndim,nel,nopen
      real*8 et
      real*8 eps(ndim)
      real*8 g

      real*8 focca(ndim),foccb(ndim)
      real*8 efa,efb
      real*8 ga,gb   
      real*8 nfoda,nfodb
      integer ihomoa,ihomob

!     convert restricted occ first to alpha/beta             
      if(nel.gt.0) then
         call occu(ndim,nel,nopen,ihomoa,ihomob,focca,foccb)
      else
         focca=0.0d0
         foccb=0.0d0
         ihomoa=0
         ihomob=0
      endif
      if(ihomoa+1.le.ndim) then 
         call FERMISMEAR(.false.,ndim,ihomoa,et,eps,focca,nfoda,efa,ga)
      endif
      if(ihomob+1.le.ndim.and.nel.gt.1) then
         call FERMISMEAR(.false.,ndim,ihomob,et,eps,foccb,nfodb,efb,gb)
      endif

      g = ga + gb

end      
