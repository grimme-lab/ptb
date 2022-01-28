!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

program gTB    
   
      use bascom   ! AO basis 
      use cbascom  ! core AO basis
      use parcom   ! TB method parameters
      use com      ! general stuff         
      use mocom    ! ref MOs for fit and momatch value

      use iso_fortran_env, only : wp => real64
      implicit none

      real(wp),allocatable :: xyz(:,:),rab(:),z(:), wbo(:,:), cn(:)
      real(wp),allocatable :: psh(:,:),q(:), psh_ref(:,:), q_ref(:), wbo_ref(:,:), qd4(:)
      real(wp),allocatable :: psh_gtb(:,:),q_gtb(:), wbo_gtb(:,:)
      real(wp),allocatable :: S(:),T(:),P(:),SS(:),F(:),D3(:,:)
      real(wp),allocatable :: eps(:),focc(:),xnorm(:)
      real(wp),allocatable :: Pref(:)
      real(wp),allocatable :: fragchrg_eeq(:),fragchrg_ref(:),fragchrg(:)
      real(wp),allocatable :: dipgrad(:,:),dipgrad_ref(:,:)
      real(wp),allocatable :: fdgrad(:,:,:),fdgrad_ref(:,:,:)
      real*4  ,allocatable :: ML1(:,:),ML2(:,:)

      integer, allocatable :: at(:)
      integer ,allocatable ::molvec(:)

      integer n
      integer ndim
      integer nopen
      integer na,nb,nel,ihomo
      integer prop
      integer ia,ib,ish,jsh,ksh,ata,atb,ati
      integer ii,jj,ij,ll
      integer molcount,idum(100),tenmap(6)
      integer i,j,k,l,m,ns,nf,nn,lin,llao(4)
      data llao/1,3,5,7 /
      real(wp) chrg ! could be fractional for model systems

      real(wp) t0,t1,w0,w1,t00,w00,ddot
      real(wp) etot,eref,enuc,ekinref,ekin,eel
      real(wp) norm,ff,f2,f1,x,y,qi,ge,r
      real(wp) dcal,dref,deeq,dang
      real(wp) a1,a2,s8
      real(wp) erfs     
      real(wp) pnt(3),dip(3),dip_ref(3),dipr(3),dipl(3)
      real(wp) alp(6),alp_ref(6),alpr(6),alpl(6)
      real(wp) efield(3),eftmp(3,6),beta(6,3)
      real(wp) floats(10),edum(86)
      real(wp),parameter :: zero = 0_wp
      character*2 asym
      character*80 str(10)
      character*80 atmp,arg1,fname,pname
      logical ex,fail,wrapo,test,test2,exref,dgrad,fdgr,raman,betaref,energ,stda,acn,rdref,nogtb
      integer TID, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM, nproc

      call timing(t00,w00)           

      rdref   =.false.
      wrapo   =.false.
      test    =.false.
      exref   =.false.
      dgrad   =.false.
      raman   =.false.
      betaref =.false.
      energ   =.false.
      stda    =.false.
      acn     =.false.
      nogtb   =.false.
      eref = 0
      ekinref = 0
      prop = 1
      pnt  = 0
      chrg = 0
      alp_ref(1)=-99
      pname='~/.atompara'

!!$OMP PARALLEL PRIVATE(TID)
!      TID = OMP_GET_THREAD_NUM()
!      IF (TID .EQ. 0) THEN
!         nproc = OMP_GET_NUM_THREADS()
!         PRINT *, '============================='
!         PRINT *, ' # OMP threads =', nproc   
!         PRINT *, '============================='
!      END IF
!!$OMP END PARALLEL 

      expscal=0
      expscal(4,1:10,1:86)=1.0d0  ! back to standard exp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! input options
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call head
      call getarg(1,fname)
      do i=2,8
      call getarg(i,arg1)
      if(index(arg1,'-avcn' ) .ne.0)  acn=.true.  ! 
      if(index(arg1,'-apo' )  .ne.0)wrapo=.true.  ! just printout shell pop for coding
      if(index(arg1,'-polar') .ne.0)prop =2       ! polar      
      if(index(arg1,'-alpha') .ne.0)prop =2       ! polar      
      if(index(arg1,'-beta') .ne.0) prop =3       ! hyperpolar      
      if(index(arg1,'-hyperpolar') .ne.0) prop =3 ! hyperpolar      
      if(index(arg1,'-energy') .ne.0)energ=.true. ! energy          
      if(index(arg1,'-stda')   .ne.0)stda=.true.  ! stda write      
      if(index(arg1,'-test').ne.0) test =.true.   ! more data output
      if(index(arg1,'-nogtb').ne.0) nogtb =.true. ! 
      if(index(arg1,'-par').ne.0)then          
      call getarg(i+1,pname)
      endif
      if(index(arg1,'-chrg').ne.0)then          
      call getarg(i+1,atmp)
      call readline(atmp,floats,str,ns,nf)
      chrg=floats(1)
      endif
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read parameter file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      write(*,*) pname
      open(unit=1,file=pname)            
      read(1,*) glob_par (1:10)
      read(1,*) glob_par(11:20)
      do i=1,44 ! change here for new elements
         read(1,*) j
         read(1,*) expscal  (1,1:10,j)   ! 1 -10  E    
         read(1,*) ener_par1 (1:10,j)    ! 11-20  E        
         read(1,*) ener_par2 (1:10,j)    ! 21-30  E   
         read(1,*) floats    (1:10)      ! 31-40  
         read(1,*) floats    (1:10)      ! 61-70  
         read(1,*) floats    (1:10)      ! 61-70  
         read(1,*) floats    (1:10)      ! 61-70  
         read(1,*) shell_xi  (1:10,j)    ! 71-80  shellQ
         read(1,*) shell_cnf1(1:10,j)    ! 81-90    "
         read(1,*) shell_cnf2(1:10,j)    ! 91-100   "
         read(1,*) shell_cnf3(1:10,j)    ! 101-110  "
         read(1,*) expscal  (3,1:10,j)   ! 111-120  "
         read(1,*) shell_cnf4(1:10,j)    ! 121-130  "
         read(1,*) shell_resp(1:10,j,1)  ! 131-140  "
         read(1,*) shell_resp(1:10,j,2)  ! 141-150  "
!        write(*,*) j
!        write(*,'(20F15.10)') shell_cnf1(1:5,j)*shell_cnf4(9,j),shell_cnf4(10,j)/shell_cnf4(9,j)
      enddo
      close(1)

! mol. charge
      inquire(file='.CHRG',exist=ex)
      if(ex)then
         open(unit=1,file='.CHRG')
         read(1,'(a)')atmp
         close(1)
         call readline(atmp,floats,str,ns,nf)
         chrg=floats(1)
      endif
! electric field
      inquire(file='.EFIELD',exist=ex)
      efield=0
      if(ex)then
         open(unit=1,file='.EFIELD')
         read(1,'(a)')atmp
         close(1)
         call readline(atmp,floats,str,ns,nf)
         if(nf.lt.3) stop '.EFIELD read error'
         efield(1:3)=floats(1:3) 
         write(*,'(''.EFIELD :'',3f12.6)') efield  
      endif
      nopen=0
!     inquire(file='.UHF',exist=ex)
!     if(ex)then
!        open(unit=1,file='.UHF')
!        read(1,'(a)')atmp
!        close(1)
!        call readline(atmp,floats,str,ns,nf)
!        nopen=int(floats(1))
!     endif

! how many atoms?
      call rd0(fname,n)

      allocate(at(n),xyz(3,n),z(n),q(n),cn(n),rab(n*(n+1)/2),        &
     &         psh(10,n),psh_ref(10,n),q_ref(n),qd4(n),wbo(n,n),     &
     &         wbo_ref(n,n),wbo_gtb(n,n),q_gtb(n),psh_gtb(10,n),     &
     &         dipgrad(3,3*n),dipgrad_ref(3,3*n),                    &
     &         fdgrad(3,3*n,9), fdgrad_ref(3,3*n,9) )

! read coordinates
      call rd(.true.,fname,n,xyz,at)
      call calcrab(n,at,xyz,rab)
 
      if(acn)then
         avcn = 0
         edum = 0
         call ncoord_erf(n,at,rab,-2.0d0,cn)
         do i=1,n 
            j=at(i)
            avcn(j)=avcn(j)+cn(i)
            edum(j)=edum(j)+1
         enddo
         open(unit=1,file='.data')
         do i=1,86
            write(1,*) avcn(i)/(edum(i)+1.d-6)
         enddo
         close(1)
         goto 9999
      endif
            
      call setavcn   ! av. el. CNs with erfs=-2.0

      idum=0
      do i=1,n
         z(i)=valel(at(i))
         idum(at(i))=idum(at(i))+1
      enddo
      nel=int(sum(z))-int(chrg)

      ndim=0
      call rdbas                      ! file: ~/.basis_vDZP
      write(*,*) 'basis read done.'
      call setupbas0(n,at,ndim)   

      allocate(S(ndim*(ndim+1)/2),P(ndim*(ndim+1)/2), SS(ndim*(ndim+1)/2), &
     &         T(ndim*(ndim+1)/2),F(ndim*(ndim+1)/2), D3(ndim*(ndim+1)/2,3), &
     &         ML1(ndim,ndim),ML2(ndim,ndim), xnorm(ndim),focc(ndim),eps(ndim),&
     &         epsref(ndim),Pref(ndim*(ndim+1)/2),cmo_ref(ndim,ndim))  

! valence basis
      call setupbas (n,at,ndim)
      write(*,*) 'basis setup done. Ndim',ndim

! core basis
      call setupcbas0(n,at)   
      call setupcbas (n,at)   

! determine occupations
      call occ(ndim,nel,nopen,ihomo,na,nb,focc)
      
! exact S and T
      call stint(n,ndim,at,xyz,rab,S,T,xnorm)    

! read DFT reference
      pnt = 0
      call rdtm(n,ndim,ihomo,at,S,focc,ekinref,dip_ref,alp_ref,Pref,rdref)
      if(rdref.and.mod(nel,2).ne.0) stop 'open-shell DFT ref. not implemented'

      if(nogtb) then
         call mlpop14(ndim,S,ML1,ML2) ! ML precalc, x=1/4
         goto 888
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! DIPGRAD
      if((.not.energ).and.(.not.stda))then
         inquire(file='dipgrad',exist=ex)
         if(ex)then
         dgrad=.true.
include 'dipgrad.f90'
         call calcrab(n,at,xyz,rab)
         call stint(n,ndim,at,xyz,rab,S,T,xnorm)    ! exact S and T
         endif
! RAMAN
         inquire(file='polgrad',exist=ex)
         if(ex)then
         raman=.true.
include 'polgrad.f90'
         call calcrab(n,at,xyz,rab)
         call stint(n,ndim,at,xyz,rab,S,T,xnorm)    ! exact S and T
         efield = 0
         endif
      endif

      if(alp_ref(1).gt.0.and.prop.lt.2) prop = 2 ! switch on polar calc
      inquire(file='beta_red',exist=ex)
      if(ex)    prop=3
      if(energ) prop=0
      if(stda)  prop=4
! SINGLE POINT
      if(prop.gt.0) call dipint(n,ndim,at,xyz,rab,xnorm,pnt,D3)! dipole integrals 
      call pgtb(.true.,prop,n,ndim,nel,nopen,ihomo,at,chrg,xyz,z,rab,pnt,xnorm,S,T,D3,&
     &          efield,qd4,ML1,ML2,psh,q,P,F,eps,eel,wbo,dip,alp) ! gd4 = EEQ to avoid recompute
      call energy(ndim,T,P,ekin) 
      write(*,'('' kinetic energy (calc/ref) :'',2f12.6)') ekin,ekinref

! ENERGY ONLY FROM P and q
      if(energ.and.(.not.stda)) then
         call gtbenergy(.true.,n,ndim,at,z,xyz,rab,q,psh,P,wbo,eref,etot)
         open(unit=12,file='.data')  
!        write(12,*) eref/float(n),etot/float(n)
         write(12,*) eref,etot
         close(12)
         open(unit=111,file='energy')
         write(111,'(''$energy'')')   
         write(111,'('' 1      '',f16.8,''   99.9 99.9 99.9 99.9'')') etot
         write(111,'(''$end'')')   
         close(111)
         goto 9999
      endif
! BETA 
!     if(prop.eq.3.or.ex)then
!        prop=3
!        if(ex)betaref=.true.
!        y=0.0005_wp ! field strength
!        write(*,'(''beta response ff-strength :'',3f10.5)') y 
!        alp = 0
!        do m=1,3
!           efield=0_wp           
!           efield(m)=y           
!           write(*,'(''Efield, dir :'',3f8.5,i4)')efield,m
!           call pgtb(.false.,102,n,ndim,nel,nopen,ihomo,at,chrg,xyz,z,rab,pnt,xnorm,S,D3,&
!    &               efield,qd4,ML1,ML2,psh,q,P,F,eps,eel,wbo,dipr,alpr) ! ALPHA calc
!           efield(m)=-y          
!           write(*,'(''Efield, dir :'',3f8.5,i4)')efield,m
!           call pgtb(.false.,102,n,ndim,nel,nopen,ihomo,at,chrg,xyz,z,rab,pnt,xnorm,S,D3,&
!    &                efield,qd4,ML1,ML2,psh,q,P,F,eps,eel,wbo,dipl,alpl) ! other direction
!           beta(1:6,m)=(alpr(1:6)-alpl(1:6))/(2_wp*y)
!           alp=alp+0.5_wp*(alpr+alpl)/3_wp ! very good approx. to true alpha (avoids extra call)
!!          call prpolar(alpr) 
!!          call prpolar(alpl) 
!        enddo
!     endif

!     compute NCI fragments (for artificial CT check)
      allocate(molvec(n))
      call mrec(molcount,xyz,n,at,molvec)
      allocate(fragchrg(molcount),fragchrg_eeq(molcount),fragchrg_ref(molcount))
      fragchrg_ref=0
      fragchrg_eeq=0
      do i=1,n
         fragchrg_ref(molvec(i))=fragchrg_ref(molvec(i))-q_ref(i)+z(i)
         fragchrg_eeq(molvec(i))=fragchrg_eeq(molvec(i))-q    (i)+z(i) 
      enddo
      if(molcount.gt.1)then
      write(*,*) 'CT ref',fragchrg_ref
      write(*,*) 'CT cal',fragchrg_eeq
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! OUTPUT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,'(/,''     ------  P-gTB results ------'')')
      write(*,'(''  # element  Z   pop.                shell populations'')')
      do i=1,n    
         ns = bas_nsh(at(i))
         floats(1:ns)=psh(1:ns,i)
         write(*,'(i3,5x,a2,f5.1,f8.4,5x,10f7.3,5x)') i,asym(at(i)),z(i),q(i),floats(1:ns)
      enddo
      call prwbo(n,at,wbo)
      call prdipole(dip)
      if(prop.eq.3) write(*,'(''polarizability computed with beta-parameterized response correction!'')')
      if(prop.eq.2) call prpolar(alp)
      if(prop.eq.3) call prbeta (beta)

 888  continue
! DFT reference exists
      if(rdref)then
      call mlpop2(n, ndim, Pref, ML1, ML2, q_ref, psh_ref)
      call wiberg(n, ndim, at, rab, Pref, S, wbo_ref)
      write(*,'(/,''     ------ reference data ------'')')
      write(*,'(''  # element  Z   pop.                shell populations'')')
      if(nogtb) open(unit=133,file='nogtb.out')
      do i=1,n    
         ns = bas_nsh(at(i))
         q_ref(i)=sum(psh_ref(1:ns,i))
         floats(1:ns)=psh_ref(1:ns,i)
         if(wrapo)call qshnorm(at(i),z(i),ns,floats)
         write(*,'(i3,5x,a2,f5.1,f8.4,5x,10f7.3)') i,asym(at(i)),z(i),q_ref(i),floats(1:ns)
         if(nogtb.and.i.eq.1) write(133,'(2i3,5x,a2,f8.4,5x,10f7.3)') i,at(i),asym(at(i)),z(i)-q_ref(i),floats(1:ns)
      enddo
      if(wrapo) stop

      call prwbo(n,at,wbo_ref)
      call prdipole(dip_ref) 
      if(alp_ref(1).gt.0) call prpolar(alp_ref) 
      if(nogtb) goto 9999

      if(test) then
         open(unit=112,file='.datap')  
         open(unit=113,file='.dataq')  
         open(unit=114,file='.datab')  
         open(unit=115,file='.datad')  
         do i=1,n
         write(113,*) z(i)-q_ref(i),z(i)-q(i) 
         do j=1,bas_nsh(at(i))
            write(112,*) psh_ref(j,i),psh(j,i) 
         enddo
         enddo
         close(112)
         close(113)
         do i=2,n
         do j=1,i-1 
            if(abs(wbo_ref(j,i)).gt.0.1) write(114,'(2F16.8,4i4)') wbo_ref(j,i), wbo(j,i), j, i, at(j),at(i)
         enddo
         enddo
         close(114)
         do i=1,3
         if(abs(dip_ref(i)).gt.0.01) write(115,*) dip_ref(i),dip(i)
         enddo
         close(115)
         if(raman)then
         open(unit=126,file='.datara')  
         do m=1,6
         do i=1,n
         do j=1,3
         if(abs(fdgrad_ref(j,i,m)).gt.1d-5) write(126,'(2F22.14)') fdgrad_ref(j,i,m), fdgrad(j,i,m)
         enddo
         enddo
         enddo
         close(126)
         endif
         if(dgrad)then
         open(unit=116,file='.datadg')  
         do i=1,3*n
         do j=1,3
         if(abs(dipgrad_ref(j,i)).gt.1d-5) write(116,'(2F22.14)') dipgrad_ref(j,i), dipgrad(j,i)
         enddo
         enddo
         close(116)
         open(unit=116,file='.datadg2')  
         do i=1,3*n
         do j=1,3
         write(116,'(2F22.14)') dipgrad_ref(j,i)
         enddo
         enddo
         write(116,*)
         do i=1,3*n
         do j=1,3
         write(116,'(2F22.14)') dipgrad(j,i)
         enddo
         enddo
         close(116)
         if(alp_ref(1).gt.0)then
         open(unit=117,file='.datapol')  
         do i=1,6  
         if(abs(alp_ref(i)).gt.1d-2) write(117,'(2F22.14)') alp_ref(i), alp(i)
         enddo
         close(117)
         endif
         endif
         if(betaref)then
         open(unit=119,file='beta_red')  
         open(unit=118,file='.databet')  
         open(unit=120,file='.databet2')  
         read(119,'(a)') atmp 
         do i=1,3
            do j=1,6
               read(119,*) x
               if(abs(x).gt.0.1) write(118,*) x,beta(j,i)
                                 write(120,*) x
            enddo
         enddo
         write(120,*)
         do i=1,3
            do j=1,6
                                 write(120,*) beta(j,i)
            enddo
         enddo
         close(118)
         close(119)
         close(120)
         endif
         goto 9999
      endif

!     call system('pwd > tmpgtb')
!     open(unit=43,file='tmpgtb')
!     read(43,'(a)') atmp 
!     close(43,status='delete')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! important fit data output file write section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      open(unit=12,file='.data')  
      if((.not.betaref).and.(.not.stda))then
! Ekin
      if(abs(ekinref) .gt.1d-8) write(12,'(2F28.14,1x,a)') 0.3d0*ekinref,0.3d0*ekin    
! momatch
      if(abs(totmatch).gt.1d-8) write(12,'(2F28.14,1x,a)') zero, 0.04d0*totmatch 
! psh
      do i=1,n
         do j=1,bas_nsh(at(i))
            write(12,'(2F28.14,1x,a)') psh_ref(j,i),psh(j,i) !,trim(atmp)
         enddo
      enddo
! mu
      f2=0.20 ! 0.3 = 15 % HCNO
      do i=1,3
         write(12,'(2F28.14,1x,a)') f2*dip_ref(i),f2*dip(i) !,trim(atmp)
      enddo
! wbo
      f2=0.25 ! 0.25=10 % HCNO
      do i=2,n
         do j=1,i-1 
            if(abs(wbo_ref(j,i)).gt.0.01) write(12,'(2F28.14,1x,a)') wbo_ref(j,i)*f2, wbo(j,i)*f2 !,trim(atmp)
         enddo
      enddo
! alpha grad
      if(raman)then
      f2=0.02       ! 0.02
      do m=1,6
      do i=1,n
         do j=1,3
            if(abs(fdgrad_ref(j,i,m)).gt.1d-3) write(12,'(2F28.14,1x,a)') fdgrad_ref(j,i,m)*f2, fdgrad(j,i,m)*f2 !,trim(atmp)
         enddo
      enddo
      enddo
      endif
! dipole grad
      if(dgrad)then
      k=0
      do i=1,n
         f2=0.65     ! 0.65 = 30 % HCNO
         if(at(i).eq.1) f2=1.10
         do ii=1,3
            k=k+1
            do j=1,3
            if(abs(dipgrad_ref(j,k)).gt.1d-4) write(12,'(2F28.14,1x,a)') dipgrad_ref(j,k)*f2, dipgrad(j,k)*f2 !,trim(atmp)
            enddo
         enddo
      enddo
      endif
! alpha
      if(alp_ref(1).gt.0)then
         f2=0.008 ! 0.01
         do i=1,6  
         if(abs(alp_ref(i)).gt.1d-2) write(12,'(2F28.14)') alp_ref(i)*f2, alp(i)*f2
         enddo
      endif
! beta
      else 
!        f2=0.1
!        open(unit=119,file='beta_red')  
!        read(119,'(a)') atmp 
!        do i=1,3
!           do j=1,6
!              read(119,*) x
!              y = beta(j,i)
!              write(12,'(2F28.14)') f2*sign(sqrt(abs(x)),x), f2*sign(sqrt(abs(y)),y)
!           enddo
!        enddo
!        close(119)
         write(*,*) 'sTDA fit data'
!        do i=1,fitcount
!           write(12,'(2F28.14)') fitdat(i)
!        enddo
      endif

      close(12)
      goto 9999
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! .data file written
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     call molecho(6,n,nel,na,nb,xyz,rab,at,z)

 9999 continue

                                          call timing(t1,w1)           
                                          call prtime(6,t1-t00,w1-w00,'all')

end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine head
      implicit none
      character(len=40),parameter:: date='Wed Jan 19 09:30:34 CET 2022'
      character(len=10),parameter:: version='3.5'

      write(*,*)
      write(*,'(7x,''=============================================='')')
      write(*,'(7x,''|                 g T B                      |'')')
      write(*,'(7x,''|                S.Grimme                    |'')')
      write(*,'(7x,''|          Universitaet Bonn, MCTC           |'')')
      write(*,'(7x,''=============================================='')')
      write(*,'(7x,''Version '',a,'', '',a)')trim(version),trim(date)
      write(*,*)
      write(*,'(7x,''Cite work conducted with this code as'')')
      write(*,'(7x,''S. Grimme, M. Mueller, A. Hansen, unpublished.'')')

      write(*,*)

end
