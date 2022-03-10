subroutine rdtm(n,ndim,homo,at,S,focc,ekin,etot,dip,alp,P,rdref)
   use iso_fortran_env, only : wp => real64
   use mocom ! fit only
   implicit none

!! ------------------------------------------------------------------------
!  Input
!! ------------------------------------------------------------------------
   integer, intent(in)    :: n                  ! number of atoms
   integer, intent(in)    :: ndim               ! number of SAOs       
   integer, intent(in)    :: homo               ! as the name says       
   integer, intent(in)    :: at(n)              ! element numbers
   real(wp),intent(in)    :: S(ndim*(ndim+1)/2) ! exact overlap maxtrix in SAO
   real(wp),intent(in)    :: focc(ndim)         ! occupations
!! ------------------------------------------------------------------------
!  Output
!! ------------------------------------------------------------------------
   real(wp), intent(out)  :: ekin
   real(wp), intent(out)  :: etot
   real(wp), intent(out)  :: dip(3)
   real(wp), intent(out)  :: alp(6)
   real(wp), intent(out)  :: P(ndim*(ndim+1)/2) ! exact P maxtrix in SAO
   logical , intent(out)  :: rdref

!! ------------------------------------------------------------------------
!  local
!! ------------------------------------------------------------------------
   real(wp)             :: xx(10), x
   real(wp),allocatable :: stmp(:,:)
   integer              :: nn, iret
   integer              :: lin        
   integer              :: i, j, k 
   character*80         :: atmp
   logical              :: pr
   logical              :: ex

   pr    = .false.
   rdref = .false.

   inquire(file='energy', exist=ex)
   if(.not.ex) return
   open(unit=10,file='energy')
10 read(10,'(a)',end=20) atmp
   call readl(atmp,xx,nn)
   if(nn.ge.2) ekin=xx(3)
   goto 10
20 continue
   close (10) 

   inquire(file='control', exist=ex)
   if(.not.ex) return
   open(unit=10,file='control')
25 read(10,'(a)',end=30) atmp
   if(index(atmp,'subenergy').ne.0) then
      read(10,'(a)',end=30) atmp
      call readl(atmp,xx,nn)
      if(nn.ge.3) etot=xx(1)
   endif
   if(index(atmp,'$dipole').ne.0) then
      read(10,'(a)',end=30) atmp
      call readl(atmp,xx,nn)
      if(nn.ge.3) dip(1:3)=xx(1:3)
   endif
   if(index(atmp,'electronic polarizability  1/3').ne.0) then
      read(10,'(a)',end=30) atmp
      call readl(atmp,xx,nn)
      alp(1)=xx(1)
      read(10,'(a)',end=30) atmp
      call readl(atmp,xx,nn)
      alp(2:3)=xx(1:2)
      read(10,'(a)',end=30) atmp
      call readl(atmp,xx,nn)
      alp(4:6)=xx(1:3)
   endif
   goto 25
30 continue
   close (10) 

   inquire(file='mos', exist=ex)
   if(.not.ex) return
   open(unit=10,file='mos',iostat=iret)
   if(iret.ne.0) return
   read(10,'(a)',end=40) atmp
   read(10,'(a)',end=40) atmp
   read(10,'(a)',end=40) atmp
   i = 0
35 read(10,'(a)',end=40) atmp
   if(index(atmp,'$end').ne.0) goto 40
   i = i + 1
   call readl(atmp,xx,nn)
   epsref (i) = xx(2)
   read(10,'(4d20.14)') cmo_ref(1:ndim,i)
   goto 35
40 continue
   close (10) 
   if(i.ne.ndim) stop 'mo counter .ne. ndim'

   allocate(stmp(ndim,ndim))
!  call prmat(6,cmo_ref,ndim,ndim,'C')
   call reordertm(ndim,n,at,stmp)

!  check proper orthonormalization for first 5 MOs
   x = 0
   do i=1,min(5,ndim)
      do j=1, ndim
         do k = 1, ndim
            x = x + cmo_ref(j,i)*S(lin(k,j))*cmo_ref(k,i)
         enddo
      enddo
   enddo
   if(abs(x-5d0).gt.1d-6) stop 'TM mos are bogous'

   call dmat(ndim,focc,cmo_ref,stmp)
   call packsym(ndim,stmp,P)

   if (pr) then
      write(*,*) 'norm :', x/5d0
      write(*,*) 'Ekin :', ekin
      write(*,*) 'dip  :', dip 
      write(*,*) 'alp  :', alp 
      write(*,*) 'eps  :', epsref
   endif

   rdref = .true.

   end

!! ------------------------------------------------------------------------
!  reorder TM mos
!! ------------------------------------------------------------------------

subroutine reordertm(nao,nat,at,stmp)
      use bascom
      use mocom  
      implicit none          
      integer, intent(in)  :: nao,nat,at(nat)

      integer i,j,ii,jj,ish,ijij,ij,lin
      integer isao,iat,ishtyp,mli,iperm(nao),llao2(0:3)
      data llao2/1,3,5,7 /
      real*8 stmp(nao,nao)

! permutation array 
! order in TM (define, infsao option):
! 1 -4 d0  = (-xx-yy+2zz)/sqrt(12) 
! 2 -3 d1a = xz    
! 3 -2 d1b = yz    
! 4 -1 d2a = xy    
! 5    d2b = (xx-yy)/2  
      isao=0
      do iat=1,nat
         do ish=1,bas_nsh(at(iat))
            ishtyp=bas_lsh(ish,at(iat))
            do mli=1,llao2(ishtyp)
               isao=isao+1
               iperm(isao)=isao
            enddo
            if(llao2(ishtyp).eq.5)then
                  iperm(isao-3 )=isao-4 ! dz2   2
                  iperm(isao-4) =isao   ! dx2y2 1
                  iperm(isao-2) =isao-1 ! dxy   3
                  iperm(isao-1) =isao-3 ! dxz   4
                  iperm(isao  ) =isao-2 ! dyz   5
            endif
         enddo
      enddo

      do i=1,nao
         ii=iperm(i)
         stmp(i,1:nao)=cmo_ref(ii,1:nao)
      enddo

      cmo_ref = stmp
 
!     call prmat(6,cmo_ref,nao,nao,'C')
!     stop 'reorder2'

end 

!! ------------------------------------------------------------------------
!  write TM mos
!! ------------------------------------------------------------------------

subroutine wr_tm_mos(nao,nat,nel,at,nopen,homo)
      use bascom
      implicit none          
      integer, intent(in)  :: nao,nel,nat,at(nat),nopen,homo

      integer i,j,ii,jj,ish,ijij,ij,lin,ihomoa,ihomob
      integer isao,iat,ishtyp,mli,iperm(nao),llao2(0:3)
      data llao2/1,3,5,7 /
      real*8,allocatable :: stmp(:,:), C(:,:), eps(:)
      real*8 focca(nao),foccb(nao),gapa,gapb,shift

      allocate(stmp(nao,nao),C(nao,nao),eps(nao))

! permutation array 
! order in TM (define, infsao option):
! 1 -4 d0  = (-xx-yy+2zz)/sqrt(12) 
! 2 -3 d1a = xz    
! 3 -2 d1b = yz    
! 4 -1 d2a = xy    
! 5    d2b = (xx-yy)/2  
      isao=0
      do iat=1,nat
         do ish=1,bas_nsh(at(iat))
            ishtyp=bas_lsh(ish,at(iat))
            do mli=1,llao2(ishtyp)
               isao=isao+1
               iperm(isao)=isao
            enddo
            if(llao2(ishtyp).eq.5)then
                  iperm(isao-3 )=isao-4 ! dz2   2
                  iperm(isao-4) =isao   ! dx2y2 1
                  iperm(isao-2) =isao-1 ! dxy   3
                  iperm(isao-1) =isao-3 ! dxz   4
                  iperm(isao  ) =isao-2 ! dyz   5
            endif
         enddo
      enddo
!     open(unit=42,file='gtb_tmpmos',form='unformatted')
      rewind 42
      read(42) C
      read(42) eps
!     close(42)

      do i=1,nao
         ii=iperm(i)
         stmp(ii,1:nao)=C(i,1:nao)
      enddo

      if(nopen.eq.0)then
      open(unit=68,file='mos.tmp')
      write(68,'(''$scfmo    scfconv=6   format(4d20.14)'')')
      write(68,'(''# SCF total energy is    -9999.9999999999 a.u.'')')
      write(68,'(''#'')')
      do i=1,homo
         write(68,'(i6,2x,''a      eigenvalue='',d20.14,''   nsaos='',i0)') i,eps(i),nao
         write(68,'(4d20.14)') stmp(1:nao,i)
      enddo
      write(68,'(''$end'')')
      close(68)
      else
      call occu(nao,nel,nopen,ihomoa,ihomob,focca,foccb) ! alpha/beta occupation
      gapa = eps(ihomoa+1)-eps(ihomoa)
      gapb = eps(ihomob+1)-eps(ihomob)
      shift= gapa - gapb
      open(unit=68,file='alpha.tmp')
      write(68,'(''$uhfmo_alpha   scfconv=6   format(4d20.14)'')')
      write(68,'(''# SCF total energy is    -9999.9999999999 a.u.'')')
      write(68,'(''#'')')
      do i=1,homo
         write(68,'(i6,2x,''a      eigenvalue='',d20.14,''   nsaos='',i0)') i,eps(i),nao
         write(68,'(4d20.14)') stmp(1:nao,i)
      enddo
      write(68,'(''$end'')')
      close(68)
!     do i=1,nao
!        write(*,'(i3,2f10.4)') i,focca(i),eps(i)
!     enddo
      open(unit=68,file='beta.tmp')
      write(68,'(''$uhfmo_beta   scfconv=6   format(4d20.14)'')')
      write(68,'(''# SCF total energy is    -9999.9999999999 a.u.'')')
      write(68,'(''#'')')
      eps(ihomob+1:nao) = eps(ihomob+1:nao) + shift  ! shift beta epsilon such that gap is the same as in alpha space
      do i=1,homo
         write(68,'(i6,2x,''a      eigenvalue='',d20.14,''   nsaos='',i0)') i,eps(i),nao
         write(68,'(4d20.14)') stmp(1:nao,i)
      enddo
      write(68,'(''$end'')')
      close(68)
!     do i=1,nao
!        write(*,'(i3,2f10.4)') i,foccb(i),eps(i)
!     enddo
      endif

end 

