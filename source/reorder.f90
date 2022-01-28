subroutine reorder(nao,nat,at,s)
      use bascom
      implicit none          
      integer, intent(in)  :: nao,nat,at(nat)
      real*8, intent(out)  :: s(nao*(nao+1)/2)

      integer i,j,ii,jj,ish,ijij,ij,lin
      integer isao,iat,ishtyp,mli,iperm(nao),llao2(0:3)
      data llao2/1,3,5,7 /
      real*8 stmp(nao*(nao+1)/2)

! permutation array for ORCA ordering
      isao=0
      do iat=1,nat
         do ish=1,bas_nsh(at(iat))
            ishtyp=bas_lsh(ish,at(iat))
            do mli=1,llao2(ishtyp)
               isao=isao+1
               iperm(isao)=isao
            enddo
            if(llao2(ishtyp).eq.3)then
                  iperm(isao-2)=isao   !z
                  iperm(isao  )=isao-1 !x
                  iperm(isao-1)=isao-2 !y
            endif
            if(llao2(ishtyp).eq.5)then
                  iperm(isao-2 )=isao  ! dz2
                  iperm(isao-3)=isao-1 ! dxz
                  iperm(isao-4)=isao-3 ! dyz
                  iperm(isao-1)=isao-4 ! dx2y2
                  iperm(isao  )=isao-2 ! dxy
            endif
         enddo
      enddo

      ij=0
      do i=1,nao
         ii=iperm(i)
         do j=1,i  
            ij = ij + 1
            jj=iperm(j)
            ijij=lin(jj,ii)
            stmp(ijij)=s(ij)
         enddo
      enddo

      s = stmp
 
end 

subroutine reorder2(nao,nat,at)
      use bascom
      use mocom  
      implicit none          
      integer, intent(in)  :: nao,nat,at(nat)

      integer i,j,ii,jj,ish,ijij,ij,lin
      integer isao,iat,ishtyp,mli,iperm(nao),llao2(0:3)
      data llao2/1,3,5,7 /
      real*8 stmp(nao,nao)

! permutation array for ORCA ordering
      isao=0
      do iat=1,nat
         do ish=1,bas_nsh(at(iat))
            ishtyp=bas_lsh(ish,at(iat))
            do mli=1,llao2(ishtyp)
               isao=isao+1
               iperm(isao)=isao
            enddo
            if(llao2(ishtyp).eq.3)then
                  iperm(isao-2)=isao   !z
                  iperm(isao  )=isao-1 !x
                  iperm(isao-1)=isao-2 !y
            endif
            if(llao2(ishtyp).eq.5)then
                  iperm(isao-2 )=isao  ! dz2
                  iperm(isao-3)=isao-1 ! dxz
                  iperm(isao-4)=isao-3 ! dyz
                  iperm(isao-1)=isao-4 ! dx2y2
                  iperm(isao  )=isao-2 ! dxy
            endif
         enddo
      enddo

      ij=0
      do i=1,nao
         ii=iperm(i)
         stmp(ii,1:nao)=cmo_ref(i,1:nao)
      enddo

      cmo_ref = stmp
 
!     call prmat(6,cmo_ref,nao,nao,'C')
!     stop 'reorder2'

end 

