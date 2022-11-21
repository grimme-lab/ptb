!--------------------------------------------------------
! compute total energy bei TB2 type expression and PTB
! input, SG 6/22
!        read(1,*) ener_par1 (1:10,j)    ! 1 -10       
!        read(1,*) ener_par2 (1:10,j)    ! 11-20       ish
!        read(1,*) expscal  (1,1:10,j)   ! 21-30       ish 6-10 also ish!
!        read(1,*) ener_par6 (1:10,j)    ! 31-40       ish
!        read(1,*) ener_par4 (1:10,j)    ! 41-50       ish
!        read(1,*) ener_par5 (1:10,j)    ! 51-60       
! fit more or less completed for HCNOF,Si-Cl Nov 2022
!--------------------------------------------------------

subroutine ptb_energy(pr,n,nao,nopen,at,z,xyz,rab,pa,psh,wbo,P,eref,etot)
      use parcom
      use com
      use bascom
      use cbascom, only: cnsao
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
      real*8, intent(inout)  :: wbo(n,n)                  ! original WBO (not used)
      real*8, intent(inout)  :: P(nao*(nao+1)/2)          ! PTB density
      real*8, intent(out)    :: eref,etot                 ! total energy 

      integer i,j,k,l,m,ij,kl,ati,atj,ia,ib,ish,jsh,ii,jj,lin,li,lj
      integer llao(4)
      data llao / 1,3,5,7 /
      integer, parameter :: idx(3, 3) = reshape([1, 2, 4, 2, 3, 5, 4, 5, 6], [3, 3])
      real*8 ea,erep,eel,edisp,ecoul,exc,ewbo,eq3,eaes,eac,ecp,ge
      real*8 r0i,r0j,r0ab,damp,fac,wbox
      real*8 zi,zj,ff,xk,r,t8,t9,arg,rcovij,ten(3)
      real*8 hi,hj,pol,tmp,tmp2,tmp3,keav,ssh,nel
      real*8 s8,a1,a2,s9
      real*8 t0,t1,w0,w1
      real*8 xx(10),par(5),sshx(86)
      character*255 atmp
      logical ex
      real*8,allocatable :: SS(:),Scv(:,:),Hdiag(:),xnorm(:),scal(:,:),spsh(:,:),q(:)
      real*8,allocatable :: sbo(:),cn(:),cn1(:),cnh0(:),h0fac(:),eps(:),wbo2(:,:,:)


      allocate(Hdiag(nao),cn1(n),cn(n),sbo(n),q(n),cnh0(n),h0fac(n), &
     &         SS(nao*(nao+1)/2),xnorm(nao),scal(10,nsh),spsh(10,n), &
     &         qm(n),dipm(3,n),qp(6,n),eps(nao),Scv(cnsao,nao),      &
     &         wbo2(n,n,4))                                             

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
      read(11) Scv  
      read(11) xnorm
      if(nopen.ne.0) read(11) spsh
      close(11)
 
!!!!!!!!!!!!!!
! dispersion
!!!!!!!!!!!!!!

      q = z - pa ! q from PTB pop  

!     GLOB_PAR            ! wB97X-3c values
      s8 =-0.5d0          ! 0 
      a1 = 0.7d0          ! 0.2464
      a2 = 4.85d0         ! 4.4737
      s9 = 0              ! 1
!                                  PTB q in D4       beta_1/2 (D4 orig 3.0,2.0)
      call dftd4_dispersion(at,xyz,q,1d0,s8,s9,a1,a2,6.0d0,3.7d0,edisp)  ! GLOB_PAR

      if(ex) eref = eref + edisp ! take same disp for fit (i.e. not considered in fit but in energy output)
                                 ! the fit is thus to wB97X/vDZP for energies

!!!!!!!!!!!!!!
! CNs         
!!!!!!!!!!!!!!

      call ncoord_erf(n,at,rab,-2d0,cn1) ! all erf damping argument are -2 here (-1.5 in camm.f90)

!     special "no H" CN used in Hdiag and atomic shift
      cn = 0
      do i = 2, n
      do j = 1, i-1 
         r = rab(lin(i,j))
         rcovij=abs(ener_par6(8,at(i)))+abs(ener_par6(8,at(j))) ! now not PTB but fitted values
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
         sbo(i) = tmp2  ! only used for output
      enddo

      ea = 0 
      eq3= 0 
      eac= 0
      tmp= 0
      tmp2=0
      do i=1,n
         ati = at(i)
         if(abs(ener_par1(2,ati)).lt.1.d-6) stop 'parameter missing'
         ea = ea +   ener_par1(1,ati)*z(i)*(1d0+ener_par1(4,ati)*cn1(i)+ener_par1(9,ati)*cn(i)) ! neutral atom energy shift

         eq3 = eq3 + ener_par1(2,ati)*q(i)**3   &               ! third order diag as in GFN2
     &             + ener_par1(3,ati)*q(i)**4   &               ! new
     &             + 0.0025d0*en(ati)*q(i)**5                   ! new, only matters for charges > 1 i.e. C^2+ or O^2-, GLOB_PAR
         t8 = 0
         do k = 1,3
            do l = 1,3
            kl = idx(l,k)
            t8 = t8 + qp(kl,i)*qp(kl,i)
            enddo
         enddo
         t9 =sum(dipm(1:3,i)**2)
         eac = eac + ener_par5(9,ati)*sqrt(t9) + ener_par5(10,ati)*sqrt(t8) ! AXC as in GFN2 but NOT as square

!        tmp=0
!        do ish=1,bas_nsh(ati)
!           tmp = tmp + spsh(ish,i)
!        enddo
!        esp = esp + ener_par5(7,ati)*tmp    ! spin density correction

      enddo

!!!!!!!!!!!!!!
! Coulomb        
!!!!!!!!!!!!!!

      call shell_es(n,at,rab,q,psh,ecoul)  ! isotropic part

      call aniso_es(n,at,rab,xyz,eaes) ! multipole part
       
!!!!!!!!!!!!!!
! ECP
!!!!!!!!!!!!!!
      
      call eecp(n,nao,at,xyz,rab,xnorm,Scv,P,ecp)

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
         fac = 1d0
         if(ati.eq.1) fac = 2d0 ! H is different, GLOB_PAR
         do ish=1,bas_nsh(ati)
            tmp = ener_par2(ish,ati) + ener_par4(ish,ati)*cn1(i)**fac + ener_par6(9,ati)*cn(i)    ! similar (but not identical) to PTB
            nel = nel + psh(ish,i)
            do j=1,llao(bas_lsh(ish,ati)+1) ! AO loop
               ii = ii + 1
               Hdiag(ii) = tmp
            enddo
         enddo
      enddo

      do i=1,n
         fac = 0.5d0
         if(at(i).eq.1) fac = 0.1d0 ! H is different, GLOB_PAR
         h0fac(i) = 1d0+ener_par1(6,at(i))*cn1(i)**fac ! just precompute
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
         ssh= SS(ij) * (hi+hj)
         atj= at(ib)
         jsh= shell2ao(j)
         lj = bas_lsh(jsh,atj)
         pol= ((hi-hj))**2            ! denominator removed (unstable in fit)
         if(ia.ne.ib) then            ! different atoms
            keav= sqrt(ener_par5(li+1,ati) * ener_par5(lj+1,atj)) 
            tmp = ssh * keav * (1d0-pol*0.1113d0)   ! R dep. removed, GLOB_PAR
         else                         ! same atoms
            if(ish.ne.jsh) then       ! s-s', p-p', d-d' off-diagonal, li=lj because S=0 otherwise
               tmp = ssh * ener_par2(8+li,ati) * h0fac(ia) * (1d0+pol*0.011d0)  ! new CN dep., GLOB_PAR
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


   call shscalE3(n,at,scal)
   call modbasd(n,at,scal)       
   call sint(n,nao,at,xyz,rab,SS,xnorm) 
   call modbas(n,at,4) 

   call wiberg2(n,nao,at,rab,P,SS,wbo2)            ! AO resolved WBO with special S

!!!!!!!!!!!!!!
! repulsion
! and WBO term
!!!!!!!!!!!!!!
      sshx(1)   = 2.5d0  ! H is different
      sshx(2:86)=0.05d0

      erep=0
      ewbo=0
      k   =0
      do i=1,n
         ati=at(i)
         r0i= ener_par1(5,ati) 
         zi = z(i) * (1d0+ener_par1(7,ati)*q(i)+ener_par2(7,ati)*q(i)**2)   ! important charge correction
         do j=1,i-1
            k   = k + 1 
            atj = at(j)

            r0j = ener_par1(5,atj) 
            zj  = z(j) * (1d0+ener_par1(7,atj)*q(j)+ener_par2(7,atj)*q(j)**2)

            r0ab= sqrt(r0i * r0j)                                      ! slightly better than amean
            damp= 0.5d0*(erf(-1.44506d0*(rab(k)-r0ab)/r0ab)+1d0)       ! GLOB_PAR

            ssh = 0.5d0*(sshx(ati) + sshx(atj))
            wbox= sum(wbo2(j,i,:)) ! total BO
            t8  = abs(wbox)**ssh            
            
            tmp = 1d0 +(ener_par4(8,atj)+ener_par4(8,ati))*t8

            erep= erep + damp * zi * zj * tmp / rab(k)**2

            tmp2 = 0
            do m=1,4 ! s-s,p-p,s-p,all-d
               tmp2 = tmp2 + 0.5d0*(ener_par5(3+m,ati)+ener_par5(3+m,atj)) * wbo2(j,i,m)
            enddo

            ewbo= ewbo + tmp2 + 0.5d0*(ener_par5(8,ati)+ener_par5(8,atj)) * wbox**2
         enddo
         k = k + 1
      enddo

!!!!!!!!!!!!!!
! add up 
!!!!!!!!!!!!!!

!     ii=idint(nel+1.d-3)
!     call free_en(nao,ii,nopen,300d0,eps,ge)
!     write(*,*)   0.4d0* (1d0+glob_par(10)), (1d0-0.589*glob_par(10))

!     GLOB_PAR
      eel = eel * 0.2471477d0   ! 0.4d0 * (1d0+glob_par(10))
      ea  =  ea * 1.2250751d0   ! (1d0-0.589*glob_par(10))
      erep=erep * 0.4749181d0   ! 
      eaes=eaes * 0.485d0              
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
         write(*,*) 
         write(*,'(''  electronic H0         :'',2F16.8)')eel   
         write(*,'(''  ECP                   :'',2F16.8)')ecp  
         write(*,'(''  WBO                   :'',2F16.8)')ewbo  
         write(*,'(''  XC (non-local)        :'',2F16.8)')exc  
         write(*,'(''electronic+WBO+ECP+XC   :'',2F16.8)')ewbo+eel+ecp+exc
         write(*,'(''  ES 2nd-order          :'',2F16.8)')ecoul
         write(*,'(''  ES anisotropic        :'',2F16.8)')eaes 
         write(*,'(''  ES higher-order       :'',2F16.8)')eq3  
         write(*,'(''total ES                :'',2F16.8)')eq3 + eaes + ecoul
         write(*,'(''  atomic multipole XC   :'',2F16.8)')eac 
         write(*,'(''  atomic increments     :'',2F16.8)')ea  
         write(*,'(''total atomic            :'',2F16.8)')ea+eac
         write(*,'(''dispersion              :'',2F16.8)')edisp
         write(*,'(''nuclear repulsion       :'',2F16.8)')erep 
         write(*,*) 
         if(abs(eref).gt.1.d-9) then
         write(*,'(''total energy (calc/ref) :'',2F16.8)')etot,eref
         write(*,'(''energy error per /atom  :'',2F16.8)')(etot-eref)/float(n)
         else
         write(*,'(''total energy            :'',2F16.8)')etot
         if(nel.eq.0) etot = etot *20d0 ! try to make H+ energy zero in fit
         endif
      endif

      deallocate(qm,dipm,qp)
      end


!! ------------------------------------------------------------------------
!! ------------------------------------------------------------------------

subroutine eecp(n,nao,at,xyz,rab,norm,Scv,P,e)
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
      real*8,  intent(in)   :: Scv(cnsao,nao)      
      real*8,  intent(out)  :: e

      integer i,j,k,l,m,nl,nn,atn,jsh,llao2(0:3),ia,ib
      data llao2/1,3,5,7 /
      real*8,allocatable :: stmp(:,:), xtmp(:,:)
      real*8 ecpfac(86)

      ecpfac = 0
      ecpfac(5:86)  = 0.05958d0 ! GLOB_PAR

      e = 0

      if(cnsao.eq.0) return

      allocate(stmp(cnsao,nao),xtmp(nao,nao))

!     N^2 step
      do i=1,nao                          
         m=0 
         do nl=1,ncorelist                     ! all atoms with core
            nn=corelist(nl)
            atn=at(nn)
            do jsh=1,cbas_nsh(atn)             ! core shells of atom nn
               do l=1,llao2(cbas_lsh(jsh,atn)) ! AOs of core shell jsh
                  m = m + 1
                  stmp(m,i)= clev(jsh,atn) * Scv(m,i) * ecpfac(atn) !  ener_par1(8,atn)
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
            if(atk.eq.1) then
               dum=ener_par1(10,atk)+sqrt(cn(k))*ener_par1(8,atk) ! H is different
            else
               dum=ener_par1(10,atk)+cn(k)      *ener_par1(8,atk)
            endif
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
         scal(ish,i) = expscal(1,ish,ati) * (1d0 + ener_par4(9,ati)*(qa-0.03d0*qa**3)) ! ^3 term ensures that factor stays positive, GLOB_PAR
         if(scal(ish,i).lt.0.1)  scal(ish,i)=0.1
         if(scal(ish,i).gt.10.0) scal(ish,i)=10.0
      enddo
   enddo

end

subroutine shscalE2(n,at,scal)
   use iso_fortran_env, only : wp => real64
   use parcom
   use bascom
   implicit none 
   integer, intent(in)    :: n                  ! number of atoms 
   integer, intent(in)    :: at(n)              ! ordinal number of atoms
   real(wp),intent(out)   ::scal(10,n)          ! exponent scaling factors
   integer  :: i,ish,ati

   do i=1,n
      ati = at(i)
      do ish=1,bas_nsh(ati)
         scal(ish,i) = expscal(1,ish+5,ati) ! * (1d0-0.005*qshtmp) ! has to be moved because for TM  elements not working!
      enddo
   enddo

end

subroutine shscalE3(n,at,scal)
   use iso_fortran_env, only : wp => real64
   use parcom
   use bascom
   implicit none 
   integer, intent(in)    :: n                  ! number of atoms 
   integer, intent(in)    :: at(n)              ! ordinal number of atoms
   real(wp),intent(out)   ::scal(10,n)          ! exponent scaling factors
   integer  :: i,ish,ati

   do i=1,n
      ati = at(i)
      do ish=1,bas_nsh(ati)
         scal(ish,i) = ener_par6(10,ati) ! * (1d0+0.005*qshtmp)
      enddo
   enddo

end

!! ------------------------------------------------------------------------
!  second order term J matrix (other mean than in PTB)
!! ------------------------------------------------------------------------

subroutine setgab3(n,at,rab,q,gab)
   use bascom
   use parcom
   use com
      implicit none
      integer, intent(in)  :: n
      integer, intent(in)  :: at(n)
      real*8,  intent(in)  :: q(n)
      real*8,  intent(in)  :: rab(n*(n+1)/2)
      real*8,  intent(out) :: gab(nsh,nsh)

      integer i,j,k,ati,atj,ish,jsh,ii,jj,lin
      real*8 gish,gjsh,xk,r2,fi,fj,cok,cmn

      cok= 0.9d0 ! GLOB_PAR
      cmn=1.0d0-cok

      ii = 0
      do i=1, n
      ati = at(i)
      fi  = gam(ati) * (1d0+q(i)*ener_par2(10,ati)) ! higher order correction for shell gamma
      do ish=1, bas_nsh(ati)
         ii = ii + 1
         gish = abs(ener_par6(ish,ati) * fi)
         jj = 0
         do j=1,n
            k = lin(j,i)
            r2= rab(k)**2
            atj = at(j)
            fj  = gam(atj) * (1d0+q(j)*ener_par2(10,atj)) ! for atom j
            do jsh=1, bas_nsh(atj)
               jj = jj + 1
               if (jj.gt.ii) cycle
               gjsh = abs(ener_par6(jsh,atj) * fj)
               xk   = sqrt  (gish * gjsh)
               gab(jj,ii)= cok/sqrt(r2+1d0/(xk**2+1d-8)) + cmn/(rab(k)+1d0/(xk+1d-8))
!              gab(jj,ii)= 1d0/sqrt(r2+1d0/(xk**2 + 1d-8))
               gab(ii,jj)= gab(jj,ii)
            enddo
         enddo
      enddo
      enddo

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

      call setgab3(n,at,rab,q,gab) 

      iish = 0
      do i=1,n
         ati = at(i)
         call shellocc_ref(ati,atocc) ! ref. atomic pop.
         do ish=1,bas_nsh(ati)
            iish = iish + 1
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Fermi smearing
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!     write(*,*) focca+foccb
!     write(*,*) eps             
!     write(*,*) ga,gb           

      g = ga + gb

end      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! special Wiberg
! orbital pair resolved, i.e., 1=s-s, 2=p-p,
!                              3=s-p, 4=all d 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine wiberg2(n,ndim,at,rab,P,S,wbo2)
      use bascom
      use gtb_la, only : la_gemm
      implicit none
      integer n,ndim,at(n)
      real*8 rab(n*(n+1)/2)
      real*8 P(ndim*(ndim+1)/2)
      real*8 S(ndim*(ndim+1)/2)
      real*8 wbo2(n,n,4)

      real*8, allocatable ::Ptmp(:,:)
      real*8, allocatable ::si  (:,:)
      real*8, allocatable ::pi  (:,:)
      integer,allocatable ::lao (:) 
      real*8 xsum,tmp
      integer i,j,k,l,m,ll,lin
      integer llao(4)
      data llao /1,3,5,7 /
      integer aose(2,n)
      integer map(0:3,0:3)

      map(0,0)=1 !s-s
      map(1,1)=2 !p-p
      map(0,1)=3 !s-p
      map(1,0)=3
      map(2,0)=4 !x-d
      map(0,2)=4
      map(2,1)=4
      map(1,2)=4
      map(2,2)=4

      allocate(Ptmp(ndim,ndim),pi(ndim,ndim),si(ndim,ndim),lao(ndim))

      call blowsym(ndim,P,pi)
      call blowsym(ndim,S,si)
      call la_gemm('N','N',ndim,ndim,ndim,1.0d0,pi,ndim,si,ndim,0.0d0,Ptmp,ndim)

      m = 1    
      do i=1,n 
         aose(1,i)=m 
         do k=1,bas_nsh(at(i))
            do l=1,llao(bas_lsh(k,at(i))+1)
               lao(m)=bas_lsh(k,at(i))
               m = m + 1 
            enddo
         enddo
         aose(2,i)=m-1
      enddo

      wbo2=0
      do i=2,n 
         do j=1,i-1
         if(rab(lin(i,j)).lt.50.0)then
            do k=aose(1,i),aose(2,i)     ! AOs on atom i
               do m=aose(1,j),aose(2,j)  ! AOs on atom j
                  ll = map(lao(k),lao(m))
                  tmp= Ptmp(k,m)*Ptmp(m,k)
                  wbo2(j,i,ll)=wbo2(j,i,ll)+tmp
               enddo
            enddo
         endif
         wbo2(i,j,1:4)=wbo2(j,i,1:4)
         enddo
      enddo

end 

