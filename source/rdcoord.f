CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c read coordinates
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine rd(echo,fname,n,xyz,iat)
      implicit real*8 (a-h,o-z)
      dimension xyz(3,n),iat(n),xx(10)
      character*128 line
      character*80  s(10)
      character*3   a3   
      character*(*) fname
      logical echo, logicals(10)
      real*8 floats(10)
      integer :: ns,nf,nl

      CHARACTER*2 el1(98),el2(98)
      DATA el1/'h ','he',
     1 'li','be','b ','c ','n ','o ','f ','ne',
     2 'na','mg','al','si','p ','s ','cl','ar',
     3 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',
     4 'zn','ga','ge','as','se','br','kr',
     5 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',
     6 'cd','in','sn','sb','te','i ','xe',
     7 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',
     8 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',
     9 'au','hg','tl','pb','bi','po','at','rn',
     1 'fr','ra','ac','th','pa','u ','np','pu','am','cm','bk','cf'/
      DATA el2/'H ','HE',
     1 'LI','BE','B ','C ','N ','O ','F ','NE',
     2 'NA','MG','AL','SI','P ','S ','CL','AR',
     3 'K ','CA','SC','TI','V ','CR','MN','FE','CO','NI','CU',
     4 'ZN','GA','GE','AS','SE','BR','KR',
     5 'RB','SR','Y ','ZR','NB','MO','TC','RU','RH','PD','AG',
     6 'CD','IN','SN','SB','TE','I ','XE',
     7 'CS','BA','LA','CE','PR','ND','PM','SM','EU','GD','TB','DY',
     8 'HO','ER','TM','YB','LU','HF','TA','W ','RE','OS','IR','PT',
     9 'AU','HG','TL','PB','BI','PO','AT','RN',
     1 'FR','RA','AC','TH','PA','U ','NP','PU','AM','CM','BK','CF'/

      ich=142
      open(unit=ich,file=fname)
      if(echo)then
      write(*,*) '         ========================='
      write(*,*) '         reading ... ',trim(fname)
      write(*,*) '         ========================='
      endif

      if(index(fname,'xyz').ne.0)then
      read(ich,'(a)')line
      read(ich,'(a)')line
      do i=1,n
         read(ich,*) a3,xyz(1:3,i)
         call elem(a3,iat(i))
      enddo
      xyz = xyz / 0.52917726
      else
      read(ich,'(a)')line
      do i=1,n
         read(ich,'(a)')line
         call readline(line, floats, logicals, s, ns, nf, nl)
         if(nf.eq.3.and.ns.eq.1)then
            xyz(1:3,i)=floats(1:3)
            call elem(s(1),iat(i))
         endif
      enddo
      endif

      close(ich)
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c read coordinates, xyz only
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine rdxyz(fname,n,xyz)
      implicit real*8 (a-h,o-z)
      dimension xyz(3,n)
      character*128 line
      character*(*) fname

      ich=142
      open(unit=ich,file=fname)

      if(index(fname,'xyz').ne.0)then
      read(ich,'(a)')line
      read(ich,'(a)')line
      do i=1,n
         read(ich,*) a3,xyz(1:3,i)
      enddo
      xyz = zyz / 0.52917726
      else
      read(ich,'(a)')line
      do i=1,n
         read(ich,*) xyz(1:3,i)
      enddo
      endif

      close(ich)
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine rd0(fname,n)
      implicit real*8 (a-h,o-z)
      dimension xx(10)
      character*128 line
      character*(*) fname
      logical ex

      n=0
      inquire(file=fname,exist=ex)
      if(.not.ex) return
      ich=142
      open(unit=ich,file=fname)

      if(index(fname,'.xyz').ne.0)then
      read(ich,*) n
      else
 100  read(ich,'(a)',end=200)line
         if(index(line,'$redu').ne.0)goto 200
         if(index(line,'$user').ne.0)goto 200
         if(index(line,'$end' ).ne.0)goto 200
         call readl(line,xx,nn)
         if(nn.ne.3) goto 100
         n=n+1
      goto 100
 200  continue
      endif

      close(ich)
      return
      end

