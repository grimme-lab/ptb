!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

program gTB    
   
      use bascom   ! AO basis 
      use cbascom  ! core AO basis
      use parcom   ! TB method parameters
      use com      ! general stuff         
      use mocom    ! ref MOs for fit and momatch value
      use dftd4     

      use iso_fortran_env, only : wp => real64
      implicit none

      real(wp),allocatable :: xyz(:,:),rab(:),z(:), wbo(:,:), cn(:)
      real(wp),allocatable :: psh(:,:),q(:), psh_ref(:,:), q_ref(:), wbo_ref(:,:), qd4(:)
      real(wp),allocatable :: S(:),T(:),P(:),tmpmat(:),F(:),D3(:,:)
      real(wp),allocatable :: eps(:),focc(:),xnorm(:)
      real(wp),allocatable :: fragchrg_ref(:),fragchrg(:)
      real(wp),allocatable :: dipgrad(:,:),dipgrad_ref(:,:)
      real(wp),allocatable :: fdgrad(:,:,:),fdgrad_ref(:,:,:)
      real(wp),allocatable :: grad(:,:),grad_ref(:,:),dum(:,:)
      real*4  ,allocatable :: ML1(:,:),ML2(:,:)

      integer, allocatable :: at(:),mapping(:)
      integer ,allocatable ::molvec(:)
      integer ,allocatable :: dgen(:)
      integer ,allocatable :: ict(:,:)

      integer n
      integer ndim
      integer nopen
      integer na,nb,nel,ihomo
      integer prop
      integer ia,ib,ish,jsh,ksh,ata,atb,ati
      integer ii,jj,ij,ll,ngrad
      integer ntrans
      integer myunit
      integer molcount,idum(100),tenmap(6)
      integer i,j,k,l,m,ns,nf,nn,lin,llao(4)
      data llao/1,3,5,7 /
      real(wp) chrg ! could be fractional for model systems

      real(wp) t0,t1,w0,w1,t00,w00,ddot
      real(wp) etot,eref,enuc,ekinref,ekin,eve,everef,egtb
      real(wp) norm,ff,f2,f1,x,y,qi,ge,r,etew,edisp,step
      real(wp) dcal,dref,deeq,dang
      real(wp) a1,a2,s8
      real(wp) erfs,er,el
      real(wp) pnt(3),dip(3),dip_ref(3),dipr(3),dipl(3)
      real(wp) scndmom(3),scndmom_ref(3)
      real(wp) alp(6),alp_ref(6),alpr(6),alpl(6)
      real(wp) efield(3),eftmp(3,6),beta(6,3)
      real(wp) floats(10),edum(86)
      real(wp) trans(9,120)
      real(wp),parameter :: zero = 0_wp
      character*2 asym
      character*80 str(10)
      character*80 atmp,arg1,fname,pname,bname
      logical ex,fail,wrapo,test,test2,tmwr,dgrad,raman_fit,ok_ekin,energ,raman
      logical stda,acn,rdref,nogtb,ok,rpbe,fitshellq,ldum,lgrad
      logical calc_ptb_grad,d4only
      integer TID, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM, nproc

      call timing(t00,w00)           

      rdref   =.false.
      wrapo   =.false.
      test    =.false.
      tmwr    =.false.
      dgrad   =.false.
      raman_fit   =.false.
      raman   =.false.
      energ   =.false.
      stda    =.false.
      acn     =.false.
      nogtb   =.false.
      rpbe    =.false.
      lgrad   =.false.
      d4only  =.false.
      calc_ptb_grad = .false.
      eref    = 0
      ekinref = 0
      prop = 1
      pnt  = 0
      chrg = 0
      nopen= 0
      alp_ref(1)=-99
      pname='~/.atompara'
      bname="~/.basis_vDZP"

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
      do i = 1, command_argument_count()
         call getarg(i, arg1)
         select case(trim(arg1))
         case('-help')
            call help
            stop
         case('-version')
            call head
            stop
         end select
      end do

      if (command_argument_count() == 0) then
         call help
         error stop
      end if

      call head
      call getarg(1,fname)
      do i=2,8
      call getarg(i,arg1)
      if(index(arg1,'-clean' ).ne.0)calc_ptb_grad=.true.  ! 
      if(index(arg1,'-avcn' ) .ne.0)  acn=.true.  ! 
      if(index(arg1,'-apo' )  .ne.0)wrapo=.true.  ! just printout shell pop for coding
      if(index(arg1,'-polar') .ne.0)prop =2       ! polar      
      if(index(arg1,'-alpha') .ne.0)prop =2       ! polar      
      if(index(arg1,'-beta') .ne.0) prop =3       ! hyperpolar      
      if(index(arg1,'-hyperpolar') .ne.0) prop =3 ! hyperpolar      
      if(index(arg1,'-energy') .ne.0)energ=.true. ! energy          
      if(index(arg1,'-e'     ) .ne.0)energ=.true. ! energy          
      if(index(arg1,'-stda')   .ne.0)stda=.true.  ! stda write      
      if(index(arg1,'-tmwr')   .ne.0)tmwr=.true.  ! TM write      
      if(index(arg1,'-test').ne.0) test =.true.   ! more data output
      if(index(arg1,'-nogtb').ne.0) nogtb =.true. ! 
      if(index(arg1,'-raman').ne.0) raman =.true. ! 
      if(index(arg1,'-d4only').ne.0) d4only =.true. ! 
!     if(index(arg1,'-fitshellq').ne.0) then
!             fitshellq =.true.                   ! output file for test   
!             nogtb =.true. 
!     endif
      if(index(arg1,'-par').ne.0)then          
      call getarg(i+1,pname)
      endif
      if(index(arg1,'-bas').ne.0)then
      call getarg(i+1,bname)
      endif
      if(index(arg1,'-chrg').ne.0)then          
      call getarg(i+1,atmp)
      call readline(atmp,floats,str,ns,nf)
      chrg=floats(1)
      endif
      if(index(arg1,'-uhf').ne.0)then          
      call getarg(i+1,atmp)
      call readline(atmp,floats,str,ns,nf)
      nopen=floats(1)
      endif
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read parameter file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      inquire(file=pname,exist=ex)
      if (.not.ex) then
         print '(a)', "Error: Cannot find parameter file '"//trim(pname)//"'.", &
            & "Provide parameter file or specify location with -par option."
         error stop
      end if

      write(*,*) pname
      open(unit=1,file=pname)            
      read(1,*) glob_par (1:10)
      read(1,*) glob_par(11:20)
      do i=1,72 ! change here for new elements
         read(1,*) j
         read(1,*) ener_par1 (1:10,j)    ! 1 -10       
         read(1,*) ener_par2 (1:10,j)    ! 11-20       
         read(1,*) expscal  (1,1:10,j)   ! 21-30  
         read(1,*) ener_par6 (1:10,j)    ! 31-40  
         read(1,*) ener_par4 (1:10,j)    ! 41-50       
         read(1,*) ener_par5 (1:10,j)    ! 51-60       
         read(1,*) expscal  (2,1:10,j)   ! 61-70  PTB
         read(1,*) shell_xi  (1:10,j)    ! 71-80    "   
         read(1,*) shell_cnf1(1:10,j)    ! 81-90    "
         read(1,*) shell_cnf2(1:10,j)    ! 91-100   "
         read(1,*) shell_cnf3(1:10,j)    ! 101-110  "
         read(1,*) expscal  (3,1:10,j)   ! 111-120  "
         read(1,*) shell_cnf4(1:10,j)    ! 121-130  "
         read(1,*) shell_resp(1:10,j,1)  ! 131-140  "
         read(1,*) shell_resp(1:10,j,2)  ! 141-150  "
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
      inquire(file='.UHF',exist=ex)
      if(ex)then
         open(unit=1,file='.UHF')
         read(1,'(a)')atmp
         close(1)
         call readline(atmp,floats,str,ns,nf)
         nopen=int(floats(1))
      endif
      inquire(file='.RPBE',exist=ex)
      if(ex)then
         tmwr=.true.
         rpbe=.true.
      endif

! how many atoms?
      call rd0(fname,n)

      allocate(at(n),xyz(3,n),z(n),q(n),cn(n),rab(n*(n+1)/2),        &
     &         psh(10,n),psh_ref(10,n),q_ref(n),qd4(n),wbo(n,n),     &
     &         wbo_ref(n,n),grad(3,n),grad_ref(3,n),dum(3,n),        &
     &         dipgrad(3,3*n),dipgrad_ref(3,3*n),ict(n,120),dgen(n), &
     &         fdgrad(3,3*n,9),fdgrad_ref(3,3*n,9),molvec(n))


! read coordinates
      call rd(.true.,fname,n,xyz,at)
      call calcrab(n,at,xyz,rab)
      if(tmwr) call wr_control_atoms(n,at,bname) ! control atoms block for e-gtb with TM
 
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

      increase_eps_weight = .false.
      idum=0
      do i=1,n
         if (metal(at(i)).ne.0) increase_eps_weight = .true.  ! increase orbital energy weight in fit
         if (tmwr.and.shell_xi(1,at(i)).ge.0d0) then
             write(*,*) i,at(i)
             stop 'element parameter missing'
         endif
         z(i)=valel(at(i))
         idum(at(i))=idum(at(i))+1
      enddo
      nel=int(sum(z))-int(chrg)

!     for PTB-RPBE D4 only (fit)
      inquire(file='ptb_dump_0',exist=ex)
      if(ex.and.d4only) then
      open(unit=11,file='ptb_dump_0',form='unformatted')
      read(11) q
      close(11)
      qd4 = z - q  ! q from pop                s8          s9          a1          a2       beta_1/2 (orig 3.0,2.0)
      call dftd4_dispersion(at,xyz,qd4,1.0d0,glob_par(2),glob_par(1),glob_par(3),glob_par(4),glob_par(5),glob_par(6),edisp)
      write(*,'('' D4 dispersion energy      :'',f14.6)') edisp
      open(unit=124,file='.EDISP')
      write(124,'(F16.8)') edisp
      close(124)
      goto 9999
      endif

      ndim=0
      call rdbas(bname)                      ! file: ~/.basis_vDZP
      write(*,*) 'basis read done.'
      call setupbas0(n,at,ndim)   

      allocate(S(ndim*(ndim+1)/2),P(ndim*(ndim+1)/2),F(ndim*(ndim+1)/2),tmpmat(ndim*(ndim+1)/2), &
     &         D3(ndim*(ndim+1)/2,3),ML1(ndim,ndim),ML2(ndim,ndim),xnorm(ndim),focc(ndim),eps(ndim)) 

      if(n.lt.200) then
         call mrec(molcount,xyz,n,at,molvec)  ! fragments (for CT check), crashes for large systems
      else
         molcount=1
         molvec(1:n)=1
      endif
      allocate(fragchrg(molcount),fragchrg_ref(molcount))

! valence basis
      call setupbas (n,at,ndim)
      write(*,*) 'basis setup done. Ndim',ndim

! core basis
      call setupcbas0(n,at)   
      call setupcbas (n,at)   

! determine occupations
      call occ(ndim,nel,nopen,ihomo,na,nb,focc)
      write(*,*) 
      write(*,*) 'nalpha ',na  
      write(*,*) 'nbeta  ',nb  
      write(*,*) 'ntotal ',na+nb
      if(nopen.eq.0.and.na.ne.nb) nopen = na - nb ! case .UHF does not exist i.e. radical
      write(*,*) 'nopen  ',nopen
      
! exact S
      inquire(file='ptb_dump',exist=ex)
      ldum = ( .not. energ ) .or. (energ .and. (.not.ex) )      ! run it in normal case or in energy mode if dump does not exist
      if(ldum) then
       call sint (n,ndim,at,xyz,rab,S,xnorm)
       tmpmat = S
      endif
                                          call timing(t1,w1)           
                                          call prtime(6,t1-t00,w1-w00,'startup and initial S')

      if(n.eq.1) call prmat(6,S,ndim,0,'overlap matrix')

! read DFT reference
      pnt = 0
      if((.not.energ).and.(.not.stda).and.(.not.tmwr)) call rdtm(n,ndim,ihomo,at,S,focc,ekinref,eref,dip_ref,alp_ref,rdref)
      if(rdref)then
         inquire(file='.no_kinetic_energy',exist=ex)
         ok_ekin = .not. ex
      endif
!     if(rdref.and.mod(nel,2).ne.0) stop 'open-shell DFT ref. not implemented'

      if(nogtb) then
         call mlpop14(ndim,S,ML1,ML2) ! ML precalc, x=1/4
         goto 888
      endif

      if(alp_ref(1).gt.0.and.prop.lt.2) prop = 2 ! switch on polar calc
      if(energ)                         prop = 0
      if(stda)                          prop = 4
      if(tmwr)                          prop = 5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! GRADIENT
      if(prop.eq.0) then
         inquire(file='gradient',exist=ex)
         if(ex)then
         lgrad=.true.
include 'grad.f90'
         call calcrab(n,at,xyz,rab)
         S = tmpmat
         if (.not.energ) ldum=.true.
         endif
      endif
! DIPGRAD
      if(prop.gt.0.and.prop.lt.4) then
         inquire(file='dipgrad',exist=ex)
         if(ex)then
         dgrad=.true.
include 'dipgrad.f90'
         call calcrab(n,at,xyz,rab)
         S = tmpmat
         endif
! RAMAN
         inquire(file='polgrad',exist=ex)
         if(ex)then
         raman_fit=.true.
include 'polgrad.f90'
         call calcrab(n,at,xyz,rab)
         S = tmpmat
         efield = 0
         endif
      endif

      if (raman) then
         allocate(mapping(6))
         mapping(1)=1
         mapping(2)=3
         mapping(3)=6
         mapping(4)=2
         mapping(5)=4
         mapping(6)=5
         write(*,'(/,a)') "--- dALPHA/dR ---"
         x=0.005_wp
         call modbas(n,at,4) 
         do i=1,n
            do j=1,3
               xyz(j,i)=xyz(j,i)+x    
               call calcrab(n,at,xyz,rab)
               call sint(n,ndim,at,xyz,rab,S,xnorm)       ! exact S
               call dipint(n,ndim,at,xyz,rab,xnorm,pnt,D3)! dipole integrals
               call pgtb(.false.,-2,n,ndim,nel,nopen,ihomo,at,chrg,xyz,z,rab,pnt,xnorm,S,D3,&
        &               efield,ML1,ML2,psh,q,P,F,eps,wbo,dip,alpr)
               xyz(j,i)=xyz(j,i)-2_wp*x
               call calcrab(n,at,xyz,rab)
               call sint(n,ndim,at,xyz,rab,S,xnorm)       ! exact S
               call dipint(n,ndim,at,xyz,rab,xnorm,pnt,D3)! dipole integrals
               call pgtb(.false.,-2,n,ndim,nel,nopen,ihomo,at,chrg,xyz,z,rab,pnt,xnorm,S,D3,&
        &               efield,ML1,ML2,psh,q,P,F,eps,wbo,dip,alpl)
               fdgrad(j,i,1:6)=(alpr(1:6)-alpl(1:6))/(2_wp*x)
               xyz(j,i)=xyz(j,i)+x   
            enddo
            write(*,'(a,i0)') "Calculated dalpha/dr for atom ", i
         enddo
         call calcrab(n,at,xyz,rab)
         S = tmpmat
         efield = 0
         open(newunit=myunit, file="polgrad.PTB", status='REPLACE',form='FORMATTED',action='WRITE')
            write(myunit,*) "$polgrad from PTB"
            do j=1,6
                do i=1,n
                    write(myunit,*) fdgrad(1:3,i,mapping(j))
                enddo
            enddo
         close(myunit)
         write(*,'(a,/)') "written dalpha/dr in TM format to 'polgrad.PTB'"
         deallocate(mapping)
     endif


! SINGLE POINT PTB
      if(ldum) then ! run it in normal case or in energy mode if dump does not exist
       if(prop.gt.0) call dipint(n,ndim,at,xyz,rab,xnorm,pnt,D3)! dipole integrals 
       call pgtb(.true.,prop,n,ndim,nel,nopen,ihomo,at,chrg,xyz,z,rab,pnt,xnorm,S,D3,&
     &          efield,ML1,ML2,psh,q,P,F,eps,wbo,dip,alp) 
       inquire(file='ptb_dump',exist=ex)
       if (ex) call system('mv ptb_dump ptb_dump_0')  ! REQUIREMENT FOR THIS COPY PROCESS WAS NOT CLEAR
      endif

! fit case
      if(rpbe)then
         call system('egtb4') ! run TM
         goto 9999
      endif

! ENERGY ONLY 
      if(energ) then
         goto 9999   ! temporary exit for PTB-RPBE testing (ptb_energy does not work for whole PSE)
         call system('cp ptb_dump_0 ptb_dump')
         call ptb_energy(.true.,n,ndim,nopen,at,z,xyz,rab,q,psh,wbo,P,eref,etot)
         open(unit=12,file='.data')  
!        write(12,*) eref/float(n),etot/float(n)
         write(12,*) eref*10d0,etot*10d0
         if(lgrad)then
            f2 =10.0
            do i=1,n
               do j=1,3
                  write(12,*) f2*grad_ref(j,i),f2*grad(j,i)
               enddo
            enddo
         endif
         close(12)
         open(unit=111,file='energy')
         write(111,'(''$energy'')')   
         write(111,'('' 1      '',f16.8,''   99.9 99.9 99.9 99.9'')') etot 
         write(111,'(''$end'')')   
         close(111)
         goto 9999
      endif

! for RPBE-PTB
      if(prop.gt.0) then
      qd4 = z - q  ! q from pop                s8          s9          a1          a2       beta_1/2 (orig 3.0,2.0)
      call dftd4_dispersion(at,xyz,qd4,1.0d0,glob_par(2),glob_par(1),glob_par(3),glob_par(4),glob_par(5),glob_par(6),edisp)
      write(*,'('' D4 dispersion energy      :'',f14.6)') edisp
      open(unit=124,file='.EDISP')
      write(124,'(F16.8)') edisp
      close(124)
      endif

      if(rdref.and.ok_ekin)then
      call tint(n,ndim,at,xyz,rab,tmpmat,xnorm)
      call energy(ndim,tmpmat,P,ekin) 
      write(*,'('' total energy (wo disp,DFT):'',2f14.6)') eref
      write(*,'('' kinetic energy (calc/DFT) :'',2f14.6)') ekin,ekinref
      endif

      if(prop.gt.0) then
         call secint(n,ndim,at,xyz,rab,xnorm,pnt,D3)! R^2 integrals 
         call secmom(n,ndim,xyz,z,xnorm,P,D3,pnt,scndmom)
         if(rdref) call secmom(n,ndim,xyz,z,xnorm,Pref,D3,pnt,scndmom_ref)
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! OUTPUT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     compute NCI fragments (for artificial CT check)
      fragchrg=0
      do i=1,n
         fragchrg(molvec(i))=fragchrg(molvec(i))-q(i)+z(i) 
      enddo

      write(*,'(/,''     ------  PTB results ------'')')
      write(*,'(''  # element  Z   pop.                shell populations'')')
      do i=1,n    
         ns = bas_nsh(at(i))
         floats(1:ns)=psh(1:ns,i)
         write(*,'(i3,5x,a2,f5.1,f8.4,5x,10f7.3,5x)') i,asym(at(i)),z(i),q(i),floats(1:ns)
      enddo
      if(molcount.gt.1) write(*,*) 'CT cal',fragchrg
      call prwbo(n,at,wbo)
      call prdipole(dip)
      call prsec(scndmom)
      if(prop.eq.3) write(*,'(''polarizability computed with beta-parameterized response correction!'')')
      if(prop.eq.2) call prpolar(alp)
      if(prop.eq.3) call prbeta (beta)

 888  continue
! DFT reference exists
      if(rdref)then
      call mlpop2(n, ndim, Pref, ML1, ML2, q_ref, psh_ref)
      call wiberg(n, ndim, at, rab, Pref, S, wbo_ref)
!     call prmat(6,wbo_ref,n,n,'WBO DFT')
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

!     open(unit=11,file='dft_dump',form='unformatted')
!     write(11) q_ref
!     write(11) psh_ref
!     write(11) wbo_ref
!     write(11) Pref
!     close(11)

!     if(dpdq) then
!        open(unit=421,file='.datadpdq')
!        write(421,'(2i4,10F16.10)') at(1),bas_nsh(at(1)),psh_ref(1:bas_nsh(at(1)),1)  ! write for run_dpsh_dq script
!        close(421)             
!        stop
!     endif

!     if(fitshellq) then
!        open(unit=112,file='.data')  
!        q=z-q_ref ! put DFT charges into the psh model
!        call guess_qsh(n,at,z,q,psh)
!        do i=1,n
!           do j=1,bas_nsh(at(i))
!              write(112,*) psh_ref(j,i),psh(j,i)
!           enddo
!        enddo
!        close(112)
!     endif

      if(wrapo) stop

      fragchrg_ref=0
      do i=1,n
         fragchrg_ref(molvec(i))=fragchrg_ref(molvec(i))-q_ref(i)+z(i)
      enddo
      if(molcount.gt.1) write(*,*) 'CT ref',fragchrg_ref

      call prwbo(n,at,wbo_ref)
      call prdipole(dip_ref) 
      call prsec(scndmom_ref)
      if(alp_ref(1).gt.0) call prpolar(alp_ref) 

      if(nogtb) goto 9999

      if(test) then
include 'testout.f90'
      endif

!     call system('pwd > tmpgtb')
!     open(unit=43,file='tmpgtb')
!     read(43,'(a)') atmp 
!     close(43,status='delete')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! important fit data output file write section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      open(unit=12,file='.data')  
      if(rdref)then
      etew = 1.05d0
      if(ok_ekin)then ! some systems are ok but T is too bad for fit
! Ekin
        if(abs(ekinref) .gt.1d-8) write(12,'(2F28.14,1x,a)') etew*ekinref,etew*ekin !,atmp
      else
        write(*,*) 'excluding T in fit'
      endif
! momatch
      inquire(file='.no_momatch',exist=ex)
      ok = .not. ex
      if(abs(totmatch).gt.1d-8.and.ok) write(12,'(2F28.14,1x,a)') zero, 0.04d0*totmatch 
! psh
      do i=1,n
         do j=1,bas_nsh(at(i))
            write(12,'(2F28.14,1x,a)') 0.9*psh_ref(j,i),0.9*psh(j,i) !,trim(atmp)
         enddo
      enddo
! mu
      f2=0.20 ! 0.3 = 15 % HCNO
      do i=1,3
         if(abs(dip_ref(i)).gt.1d-4) write(12,'(2F28.14,1x,a)') f2*dip_ref(i),f2*dip(i) !,trim(atmp)
      enddo

! sec mom
      f2=0.02
      do i=1,3
         write(12,'(2F28.14,1x,a)') f2*scndmom_ref(i),f2*scndmom(i) !,trim(atmp)
      enddo

! wbo
      f2=0.20 ! 0.25=10 % HCNO
      x =0.01 ! cut-off
      do i=2,n
         do j=1,i-1 
            if(abs(wbo_ref(j,i)).gt.x) write(12,'(2F28.14,1x,a)') wbo_ref(j,i)*f2, wbo(j,i)*f2 !,trim(atmp)
         enddo
      enddo
! alpha grad
      if(raman_fit)then
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
         if(abs(alp_ref(i)).gt.1d-1) write(12,'(2F28.14)') alp_ref(i)*f2, alp(i)*f2
         enddo
      endif
      endif

      close(12)
      goto 9999
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! .data file written
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     call molecho(6,n,nel,na,nb,xyz,rab,at,z)

 9999 continue

      if(tmwr) call system('touch .gtbok')

                                          call timing(t1,w1)           
                                          call prtime(6,t1-t00,w1-w00,'all')

end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine help
      use, intrinsic :: iso_fortran_env, only : output_unit
      implicit none
      write(output_unit, '(a)') &
         "Usage: ptb <input> [options]...", &
         "", &
         "Geometry input must be given with the first argument.", &
         "Accepted formats are Turbomole coord and xyz format.", &
         "", &
         "Options:", &
         "", &
         "-chrg <int>        specify systems total charge", &
         "-uhf <int>         specify systems #open shells", &
         "-avcn              just calculate coordination numbers", &
         "-apo               just printout shell pop for testing", &
         "-polar/-alpha      calculate polarizibility", &
         "-hyperpolar/-beta  calculate hyperpolarizibility", &
         !"-energy            evaluate (electronic) energy", &
         "-stda              output stda/TM compatible format", &
         "-nogtb             skip gTB calculation", &
         "-par <file>        read parameters from provided file", &
         "-bas <file>        read basis set from provided file", &
         "-test              more printout for testing and debugging", &
         "-version           print version header and exit", &
         "-help              show this help message", &
         ""
end subroutine help


subroutine head
      implicit none
      character(len=40),parameter:: date='Do 28. Apr 14:28:45 CEST 2022'
      character(len=10),parameter:: version='3.7'

      write(*,*)
      write(*,'(7x,''=============================================='')')
      write(*,'(7x,''|                 P T B                      |'')')
      write(*,'(7x,''|                S.Grimme                    |'')')
      write(*,'(7x,''|          Universitaet Bonn, MCTC           |'')')
      write(*,'(7x,''=============================================='')')
      write(*,'(7x,''Version '',a,'', '',a)')trim(version),trim(date)
      write(*,*)
      write(*,'(7x,''Cite work conducted with this code as'')')
      write(*,'(7x,''S. Grimme, M. Mueller, A. Hansen, unpublished.'')')

      write(*,*)

end

