
      subroutine wr_control_atoms(n,at)  
      implicit none
      integer n, at(n)
      character*10000 z
      character*80 a80 
      CHARACTER*2 EL(107)
      DATA EL/'h ','he',
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
      integer i,j,k,l,nn,ie,ib,ierr,io,idum(86),idum2(n)
 
      idum = 0
      do i=1,n
         idum(at(i))=idum(at(i))+1
      enddo

      io=11
      open(unit=io,file='control.tmp')
      do i=1,86
         if(idum(i).gt.0) then
            idum2 = 0
            do j=1,n
               if(at(j).eq.i) idum2(j)=2
            enddo
            z=''
            call writil(n,2,1,1000,ierr,idum2,z)
            l=len(trim(z))
            nn=l / 75 + 1
            ib=1
            do l=1,nn
               ie=ib+74
  10           if(z(ie:ie).ne.','.and.l.ne.nn) then
                    ie=ie-1
                    goto 10
               endif
               if(ib.eq.1) then
                  write(a80,'(a2,1x,a,1x,''\'')') el(i), z(ib:ie)
               else
                  write(a80,'(3x   ,a,1x,''\'')')        z(ib:ie)
               endif
               ib=ie+1  
               write(io,'(a80)') a80
            enddo
            if(i.gt.1)then
            write(io,'(''   basis ='',a2,'' vDZP            \'')') el(i)
            else
            write(io,'(''   basis ='',a2,'' vDZP_small      \'')') el(i)
            endif
            if(i.gt.4.and.i.lt.11)
     &      write(io,'(''    ecp  ='',a2,'' ecp2-vDZP \'')') el(i)
            if(i.gt.10)
     &      write(io,'(''    ecp  ='',a2,'' ecp-vDZP \'')') el(i)
            write(io,'(''    jbas ='',a2,'' vDZP'')') el(i)
         endif
      enddo
      close(io)

      return
      end
       

c ---------------------------------------------------------------------
      subroutine writil(n,k,itab,ieol,ierr,idef,zeile)
      implicit real*8 (a-h,o-z)

c ---------------------------------------------------------------------
c     subroutine writes list of indices onto string <zeile> starting
c     at position <itab> of <zeile>. indices i have to obey the
c     condition idef(i) >= k. if the output would exceed beyond
c     position <ieol> of <zeile>, the return value <ierr> will be
c     set to -1 and to 0 otherwise. if it is possible, the last colon
c     will be printed. note that idef(i):=-idef(i) for all indices
c     which have already been digested.
c ---------------------------------------------------------------------

      dimension idef(*)
      character zeile*(*),zahl1*8,zahl2*8,komma,strich

      data komma,strich /',','-'/

      i=0
      ierr=-1
      ikomma=0
      ioff=itab-1
  100 i=i+1
      if(i.gt.n) goto 900
      if(idef(i).lt.k) goto 100
      istart=i
  200 i=i+1
      if(i.gt.n) goto 300
      if(idef(i).ge.k) goto 200
  300 iend=i-1
      if(istart.eq.iend) then
c     --- output of single integer
        call itostr(istart,lzahl1,zahl1)
        if(ikomma.ne.0) then
          ioff=ioff+1
          if(ioff.gt.ieol) return
          zeile(ioff:ioff)=komma
        endif
        if(ioff+lzahl1.gt.ieol) return
        zeile(ioff+1:ioff+lzahl1)=zahl1(1:lzahl1)
        ioff=ioff+lzahl1
        idef(istart)=-idef(istart)
      else
c     --- output of two integers
        call itostr(istart,lzahl1,zahl1)
        if(ikomma.ne.0) then
          ioff=ioff+1
          if(ioff.gt.ieol) return
          zeile(ioff:ioff)=komma
        endif
        call itostr(iend  ,lzahl2,zahl2)
        if(ioff+lzahl1+1+lzahl2.gt.ieol) return
        zeile(ioff+1:ioff+lzahl1)=zahl1(1:lzahl1)
        ioff=ioff+lzahl1
        zeile(ioff+1:ioff+1)=strich
        ioff=ioff+1
        zeile(ioff+1:ioff+lzahl2)=zahl2(1:lzahl2)
        ioff=ioff+lzahl2
        do 400 i=istart,iend
  400    idef(i)=-idef(i)
      endif
      ikomma=1
      goto 100
  900 ierr=0

      return
      end


c ---------------------------------------------------------------------
      subroutine itostr(n,l,a)
      implicit real*8 (a-h,o-z)

c ---------------------------------------------------------------------
c      turn integer value "n" into the string "a(1:l)"
c ---------------------------------------------------------------------

      character*(*) a
      character*10 zahlen
      data zahlen /'0123456789'/

      nnn=n
      if(n.lt.0) nnn=-n
      imax=len(a)
      iput=imax
  100 mmm=nnn/10
      mrest=nnn-mmm*10+1
      if((n.lt.0.and.iput.eq.1).or.iput.eq.0) stop '<itostr overflow>'
      a(iput:iput)=zahlen(mrest:mrest)
      iput=iput-1
      if(mmm.eq.0) goto 900
      nnn=mmm
      goto 100
  900 l=imax-iput
      if(n.lt.0) then
        a(iput:iput)='-'
        iput=iput-1
        l=l+1
      endif
      if(l.eq.imax) return
      do 200 k=1,l
  200 a(k:k)=a(iput+k:iput+k)
      do 300 k=l+1,imax
  300 a(k:k)=' '

      return
      end
