ccccccccccccccccccccccccccccccccccccccccccc
!    write out unformatted sTDA input     c
ccccccccccccccccccccccccccccccccccccccccccc           
! ncent  : # atoms
! nmo    : # MOs
! nbf    : # AOs
! nprims : # primitives (in total)
! xyz(4,ncent) : Cartesian coordinates & nuclear charge
! cont(nprims) : contraction coefficients of primitives
! alp(nprims) : exponents of primitives
! cmo(nbf,nmo) : LCAO-MO coefficients 
! eval(nmo)    : orbital eigenvalues
! occ(nmo)     : occupation # of MO
! ipty(nprims) : angular momentum of primitive function
! ipao(nbf)    : # primitives in contracted AO 
! ibf(ncent)   : # of contracted AOs on atom

      subroutine printmos(ncent,at,xyz,nmo,homo,norm,mowrcut,eval,tmp)

      use bascom
      implicit none
      
      integer, intent ( in ) :: nmo,ncent,at(ncent),homo
      real*8,  intent ( in ) :: xyz(3,ncent)
      real*8,  intent ( in ) :: norm(nmo)    
      real*8,  intent ( in ) :: mowrcut           
      real*8,  intent ( in ) :: eval(nmo)
      real*8,  intent ( inout ) :: tmp(nmo,nmo)

      ! temporary variables
      integer nbf
      integer i,j,k,nprims,nmomax,iat,iwfn,iao
      real*8 dum
      character*2 atyp
      integer lao(ncao),nprim(ncao),aoatcart(ncao)
      integer lladr(0:3),ll(0:3)
      data lladr  /1,3,6,10/
      data ll     /0,1,4,10/
      real*8,allocatable :: occ (:)
      real*8,allocatable :: cmo(:,:)

      allocate(cmo(ncao,nmo),occ(nmo))

!     rewind(42)
!     read(42) tmp
!     read(42) eval
      
      do i=1,nmo 
         tmp(i,:)=tmp(i,:)*norm(i)
      enddo

      call sao2cao(nmo,tmp,cmo,ncent,at)

      do i=1,nmo 
         tmp(i,:)=tmp(i,:)/norm(i)
      enddo

      nbf=ncao
      nprim=prim_npr
      nprims=npr

      iwfn=29
      open(unit=iwfn,file='wfn.xtb',form='unformatted',
     .     status='replace')

      occ = 0
      occ(1:homo) = 2.0d0

! only print out virtuals below cutoff
      nmomax=0
      do i=1,nmo
         if(eval(i).gt.mowrcut.and.nmomax.eq.0)nmomax=i-1
      enddo
      if(nmomax.eq.0) nmomax=nmo

                    !***********
                    ! RHF case *
                    !***********
! write dimensions
       write(iwfn)1
       write(iwfn)ncent,nbf,nmomax,nprims
! now write coordinates & atom symbol
       do i = 1,ncent
         call aasym(at(i),atyp)
         write(iwfn) atyp
       enddo

       do i = 1,ncent
         do j=1,3
            dum=xyz(j,i)
            write(iwfn) dum
         enddo 
         write(iwfn) at(i)
       enddo       
! Now print basis set data                       

      k=0
      do i=1, ncent
         iat=at(i)
         do j=1, bas_nsh(iat)
         do iao=1,lladr(bas_lsh(j,iat))
         k=k+1
         lao(k)=ll(bas_lsh(j,iat))+iao
         aoatcart(k)=i
         enddo
         enddo
      enddo

! print ipty
       do i=1,nbf
          k = lao(i)
          do j=1,nprim(i)
             write(iwfn) k
          enddo
       enddo
! iaoat
       do i=1,nbf
          k=aoatcart(i)
          do j=1,nprim(i)
             write(iwfn) k
          enddo
       enddo
! ipao
       do i=1,nbf
          k=i
          do j=1,nprim(i)
             write(iwfn) k
          enddo
       enddo       

! exponents and coefficients                                                                                                         
         write(iwfn) prim_exp(1:nprims)
         write(iwfn) prim_cnt(1:nprims) 
         
! now the mo data
         write(iwfn) occ(1:nmomax) 
!          allocate(eval_shift(nmomax))
!          Do i=1,nmomax
!          if(occ(i).lt.1.d-6)then
!          eval_shift(i)=eval(i) + 0.17d0 !shift virtuals
!          else
!          eval_shift(i)=eval(i) 
!          endif
!          enddo
!          write(iwfn) eval_shift(1:nmomax)
         write(iwfn) eval(1:nmomax)
!          deallocate(eval_shift)
  
         write(iwfn) cmo(1:nbf,1:nmomax)
      close(iwfn)

      return
      end

C     *****************************************************************         

      subroutine AASYM(I,asy)
      integer, intent ( in ) :: i
      CHARACTER*2 ASY
      CHARACTER*2 ELEMNT(107), AS
      DATA ELEMNT/'h ','he',
     1 'li','be','b ','c ','n ','o ','f ','ne',
     2 'na','mg','al','si','p ','s ','cl','ar',
     3 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',
     4 'zn','ga','ge','as','se','br','kr',
     5 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',
     6 'cd','in','sn','sb','te','i ','xe',
     7 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',
     8 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',
     9 'au','hg','tl','pb','bi','po','at','rn',
     1 'fr','ra','ac','th','pa','u ','np','pu','am','cm','bk','cf','xx',
     2 'fm','md','cb','xx','xx','xx','xx','xx'/
      AS=ELEMNT(I)
      CALL UPPER(AS)
      ASY=AS
      if(i.eq.103) asy='XX'
      RETURN
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c transforms sao(5d) integrals to cao(6d) basis 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine sao2cao(nbf,s,x,ncent,at)
      use bascom
      implicit none

      integer nbf,new,ncent,at(ncent)
      real*8  s(nbf,nbf),x(ncao,nbf)
      real*8  xcart
      integer lll(20),firstd(nbf),idprev
      integer i,j,k,jj,mm,m
      data lll/1,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4/

      real*8 trafo(5,6)
      
      integer lao(ncao),iat,iao
      integer lladr(0:3),ll(0:3)
      data lladr  /1,3,6,10/
      data ll     /0,1,4,10/

      k=0
      do i=1, ncent
         iat=at(i)
         do j=1, bas_nsh(iat)
         do iao=1,lladr(bas_lsh(j,iat))
         k=k+1
         lao(k)=ll(bas_lsh(j,iat))+iao
         enddo
         enddo
      enddo
      
      
      
! Sign of 'trafo(2,:)' changed with respect to xTB
      trafo = 0.0d0 
! x2 
      trafo(1,1)=1./dsqrt(2.d0)*dsqrt(3.d0/2.d0)
!      trafo(2,1)=0.408248290464D+00*sqrt(3./2.)
      trafo(2,1)=-0.50d0
! y2 
      trafo(1,2)=-1./dsqrt(2.d0)*dsqrt(3.d0/2.d0)
!      trafo(2,2)=0.408248290464D+00*sqrt(3./2.)
      trafo(2,2)=-0.50d0
! z2
      trafo(1,3)=0.0d0
      trafo(2,3)=1.0d0
!      trafo(2,3)=-0.816496580928D+00*sqrt(3./2.)
c rest
      trafo(3,4)=1.0d0
      trafo(4,5)=1.0d0
      trafo(5,6)=1.0d0

      new=ncao-nbf

      if(new.eq.0) then
         return
      endif

      firstd = 0
      i=1    
      j=0 
      ! lao is still in old dimensions (i.e., ncao) while s comes with nsao
 42   if(lao(i).gt.4.and.lao(i).le.10)then
         firstd(i-j:i-j+4)=i-j
         j=j+1
         i=i+5
      endif
      i=i+1
      if(i.lt.ncao)goto 42
      ! sanity check
      if(new.ne.j) stop 'error in sao2cao trafo'
     
      x=0.0d0

      do i=1,nbf ! go through eigenvectors
         k = 0
         idprev=0
         do j=1,nbf ! go through LCAO-MO coefficients
            if(idprev.gt.0.and.firstd(j).eq.idprev) cycle 
            if(firstd(j).gt.idprev)then ! if a set of d functions is found, do trafo for all six d orbitals 
               do jj=1,6
                 k=k+1
                 xcart=0.0d0
                 do m=1,5
                    mm=firstd(j)-1+m
                    xcart=xcart+trafo(m,jj)*s(mm,i)
                 enddo
                 x(k,i)=xcart
                enddo 
                idprev=firstd(j) ! setting idprev to new value guarantees that the following 4 spherical d functions will be skipped (we already did the trafo)
                cycle
            endif
            k=k+1
            x(k,i)= s(j,i)
         enddo
         if (k.ne.ncao) stop 'error in eigenvector dimension'
      enddo

      end
