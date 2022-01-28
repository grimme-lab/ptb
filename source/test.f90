subroutine calcpauli(n,nao,nsh,at,S,T,v)
      use  bascom
      use  parcom
      implicit none          
      integer, intent(in)   :: nsh,nao,n,at(n)
      real*8,  intent(out)  :: S(nao*(nao+1)/2)    
      real*8,  intent(out)  :: T(nao*(nao+1)/2)    
      real*8,  intent(out)  :: v(nao*(nao+1)/2)    

      integer i,j,k,l,m,nl,atn,jsh,llao2(0:3)
      data llao2/1,3,5,7 /
      real*8 vecp, ddot, atocc(10)
      real*8,allocatable :: stmp(:,:), sdum(nao,nao)

      allocate(stmp(nao,nao),sdum(nao,nao))
      call blowsym(nao,S,sdum)

!     N^2 step
      do i=1,nao                          
         m=0 
         do nl=1,n                           ! all atoms
            atn=at(nl)
            call shellocc_ref(atn,atocc)     ! ref. atomic pop.
            do jsh=1,bas_nsh(atn)            ! shells of atom nn
              do l=1,llao2(bas_lsh(jsh,atn)) ! AOs of shell jsh
               m = m + 1
               stmp(m,i)=-T(lin(m,m)) * S(m,i) * atocc(jsh)
              enddo
            enddo
         enddo
      enddo

!     N^3 step
      v = 0 
      k = 0 
      do i=1, nao 
         do j=1, i
            k = k + 1 
            v(k) = -ddot(nao,stmp(1,i),1,S(1,j),1)
         enddo
      enddo

      end

