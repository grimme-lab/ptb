!! ------------------------------------------------------------------------
!  MO match score for occupied space and scaled LUMO contribution
!! ------------------------------------------------------------------------
subroutine momatch(pr,ndim,nocc,nvirt,S)
      use mocom  
      use gtb_la, only : la_gemm, la_symm
      implicit none          
      logical, intent(in)  :: pr                
      integer, intent(in)  :: ndim,nocc,nvirt
      real*8 , intent(in)  :: S(ndim*(ndim+1)/2)
      
      integer i,j,k,l,match
      real*8,allocatable :: SS(:,:), C2(:,:), C(:,:), e(:), srt(:)
      real*8 norm, wei, ss2, ss4, smax, gap_ref, homo_ref
      real*8,parameter :: lumoweight = 1d0/2d0  ! LUMO weight for fit
      real*8,parameter :: mothr      = 0.3d0    ! lower weight for MOs .lt. HOMO of this vaue
      real*8,parameter :: au2ev = 27.2113957d0 

      allocate(SS(ndim,ndim),C2(ndim,ndim),C(ndim,ndim),e(ndim),srt(nocc))

      if(pr) then
         write(*,*) 'running MO match with DFT reference'
         write(*,*) ' MO  DFT #   eps (DFT)    eps      penalty       Smax(1,2)'
      endif

      rewind 42   
      read  (42) C  ! gTB MOs, DFT on cmo_ref        (see mocom)
      read  (42) e  !  "  eigenvalues, DFT on epsref   "   2

      call blowsym(ndim,S,SS)
      CALL la_symm('L','L',ndim,ndim,1.D0,SS,ndim,C,ndim,0.D0,C2,ndim)  
      call la_gemm('T','N',ndim,ndim,ndim,1.0d0,cmo_ref,ndim,C2,ndim,0.0d0,SS,ndim)

      homo_ref=epsref(nocc)
      totmatch = 0
      do i=1, nocc
         norm = 0
         smax = 0
         ss4  = 0
         match= i
         do j=1, nocc
            ss2 = SS(j,i)**2
            srt(j) = -ss2
            if(homo_ref-epsref(j).gt.mothr) then
               wei = 0.4
            else
               wei = 1.5
            endif
            norm = norm + ss2 * abs(epsref(j) - e(i)) * wei
            if(ss2.gt.smax)then ! get best matching DFT MO number for printout
               smax = ss2
               match = j
            endif
            ss4 = ss4 + ss2**2
         enddo
         if(homo_ref-epsref(i).gt.mothr) then
            norm = norm + (1d0-ss4)*0.02d0 ! add spread of overlap i.e. 10000... is best
         else
            norm = norm + (1d0-ss4)*0.08d0 ! 
         endif
         call qqsort(srt,1,nocc)
         if(pr)write(*,'(2i4,6f11.4)') i,match,epsref(match),e(i),norm,-srt(1),-srt(2)
         totmatch = totmatch + norm 
!        write(133,'(6f11.4)') epsref(match),e(i)
      enddo
!     write(133,'(6f11.4)') epsref(nocc+1),e(nocc+1)

      if(pr)write(*,*) 'MO match score occupied MOs   :', totmatch  

      gap_ref = (epsref(nocc+1) - epsref(nocc))*au2ev
      if(gap_ref.lt.10d0) then ! only reasonable gaps included in fit
         do j=nocc+1,nocc+2
          totmatch = totmatch + abs(epsref(j) - e(j)) * dble(nocc) * lumoweight  ! add LUMO deviation without assignment
         enddo                                                                   ! ie assume its right (because it should be the EA)
      endif

      if(pr)  write(*,*) 'total MO match score with gap :', totmatch  ! should be zero for perfect MOs
         
      deallocate(srt)
      allocate(srt(ndim-nocc+1))
      do i=nocc + 1, nocc + nvirt ! nvirt ca. 20
         k = 0
         ss4 = 0
         smax = 0
         do j=nocc + 1, ndim        
            ss2 = SS(j,i)**2
            k = k + 1
            srt(k) = -ss2
            ss4 = ss4 + ss2**2
            if(ss2.gt.smax)then 
               smax = ss2
               match = j
            endif
         enddo
         call qqsort(srt,1,nvirt)
         if(abs(e(i)-epsref(match)).lt.0.2d0)then
             norm = (1d0-ss4)*0.05d0 
         else
             norm = 0d0
         endif
         totmatch = totmatch + norm
         if(pr)write(*,'(2i4,6f11.4)') i,match,epsref(match),e(i),norm,-srt(1),-srt(2)
      enddo 

      if(pr)then
         write(*,*) 'total MO match score incl virt:', totmatch  
         write(*,'('' gap (eV)  DFT vs. gTB        : '',2f9.5)') gap_ref,(e(nocc+1) - e(nocc))*au2ev
      endif

      end
