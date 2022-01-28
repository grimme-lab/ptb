ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c fermi smearing      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE FERMISMEAR(prt,NORBS,NEL,T,eig,occ,fod,e_fermi,s)
      IMPLICIT NONE                        
      integer norbs 
      integer nel   
      real*8  eig(norbs)
      real*8  occ(norbs)
      real*8  t
      real*8  fod
      real*8  e_fermi
      LOGICAL PRT

      real*8 boltz,bkt,occt,total_number,thr
      real*8 total_dfermi,dfermifunct,fermifunct,s,change_fermi

      PARAMETER (BOLTZ = 3.166808578545117E-06) !*27.2113957)
      PARAMETER (thr   = 1.D-9)
      integer ncycle,i,j,m,k,i1,i2

      BKT = BOLTZ*T
 
      E_FERMI = 0.5*(EIG(NEL)+EIG(NEL+1))
      OCCT=NEL

      NCYCLE = 0
 10     TOTAL_NUMBER = 0.0
        TOTAL_DFERMI = 0.0
        NCYCLE = NCYCLE+1
        DO I = 1, NORBS
          FERMIFUNCT = 0.0
          if((EIG(I)-E_FERMI)/BKT.lt.50) then 
            FERMIFUNCT = 1.0/(EXP((EIG(I)-E_FERMI)/BKT)+1.0)
            DFERMIFUNCT = EXP((EIG(I)-E_FERMI)/BKT) /                
     .      (BKT*(EXP((EIG(I)-E_FERMI)/BKT)+1.0)**2)
          ELSE
            DFERMIFUNCT = 0.0
          END IF
          OCC(I) = FERMIFUNCT
          TOTAL_NUMBER = TOTAL_NUMBER + FERMIFUNCT
          TOTAL_DFERMI = TOTAL_DFERMI + DFERMIFUNCT
        END DO
        CHANGE_FERMI = (OCCT-TOTAL_NUMBER)/(TOTAL_DFERMI+1d-14)
        E_FERMI = E_FERMI+CHANGE_FERMI
      IF(ABS(OCCT-TOTAL_NUMBER).GT.thr.AND.NCYCLE.LT.200) GOTO 10

      fod=0
      s  =0
      do i=1,norbs
      if(occ(i).gt.thr.and.1.0D00-OCC(I).gt.thr)
     .   S=S+OCC(I)*LOG(OCC(I))+(1.0d0-OCC(I))*LOG(1.0D00-OCC(I))
         if(eig(i).lt.e_fermi)then
            fod=fod+1.0d0-occ(i)
         else
            fod=fod+      occ(i)
         endif
      enddo
      s=s*3.166808578545117E-06*t

      IF(PRT)THEN
      WRITE(*,'('' T,E(Fermi),NFOD : '',2F10.3,F10.6)') T,E_FERMI,fod
      ENDIF

      RETURN
      END

      subroutine occu(ndim,nel,nopen,ihomoa,ihomob,focca,foccb)
      implicit none
      integer nel,nopen,ndim,ihomoa,ihomob
      real*8 focca(ndim), foccb(ndim)
      integer focc(ndim)
      integer i,na,nb,ihomo

      focc=0
      focca=0
      foccb=0
c even nel      
      if(mod(nel,2).eq.0)then
      ihomo=nel/2
      do i=1,ihomo 
         focc(i)=2
      enddo
      if(2*ihomo.ne.nel) then
         ihomo=ihomo+1
         focc(ihomo)=1
         if(nopen.eq.0)nopen=1
      endif
      if(nopen.gt.1)then
         do i=1,nopen/2
            focc(ihomo-i+1)=focc(ihomo-i+1)-1
            focc(ihomo+i)=focc(ihomo+i)+1
         enddo
      endif
c odd nel      
      else
      na=nel/2+(nopen-1)/2+1
      nb=nel/2-(nopen-1)/2
      do i=1,na             
         focc(i)=focc(i)+1
      enddo
      do i=1,nb             
         focc(i)=focc(i)+1
      enddo
      endif

      do i=1,ndim
         if(focc(i).eq.2)then
            focca(i)=1.0d0
            foccb(i)=1.0d0
         endif
         if(focc(i).eq.1)focca(i)=1.0d0
      enddo

      ihomoa=0
      ihomob=0
      do i=1,ndim
         if(focca(i).gt.0.99) ihomoa=i
         if(foccb(i).gt.0.99) ihomob=i
      enddo

      if(ihomoa.lt.1) stop 'internal error in occu'
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
