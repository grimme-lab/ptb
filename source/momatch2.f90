!! ------------------------------------------------------------------------
!  MO match score for occupied space and scaled LUMO contribution
!! ------------------------------------------------------------------------
subroutine momatch2(pr,ndim,nocc,S)
      use mocom  
      use gtb_la, only : la_symm, la_gemm
      implicit none          
      logical, intent(in)  :: pr                
      integer, intent(in)  :: ndim,nocc
      real*8 , intent(in)  :: S(ndim*(ndim+1)/2)
      
      integer i,j,k,l,match
      real*8,allocatable :: SS(:,:), C2(:,:), C(:,:), e(:), srt(:)
      real*8 norm, wei, ss2, ss4, smax, gap_ref, homo_ref
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
      fitcount = 0
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
               wei = 1.2
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
            norm = norm + (1d0-ss4)*0.06d0 ! 
         endif
         fitcount = fitcount + 1
         fitdat(fitcount) = norm
         call qqsort(srt,1,nocc)
         if(pr)write(*,'(2i4,6f11.4)') i,match,epsref(match),e(i),norm,-srt(1),-srt(2)
         totmatch = totmatch + norm 
      enddo

      if(pr)write(*,*) 'MO match score occupied MOs   :', totmatch  

!     add some virt. but overlap only in penalty
      deallocate(srt)
      allocate(srt(ndim-nocc+1))
      do i=nocc + 1, ndim          
         if(epsref(i).gt.0.5d0) cycle
         fitcount = fitcount + 1
         k   = 0
         ss4 = 0
         smax= 0
         wei = 0.02
         norm= 0
         do j=nocc + 1, ndim        
            ss2 = SS(j,i)**2
            norm = norm + ss2 * abs(epsref(j) - e(i)) * wei
            k = k + 1
            srt(k) = -ss2
            ss4 = ss4 + ss2**2
            if(ss2.gt.smax)then 
               smax = ss2
               match = j
            endif
         enddo
         call qqsort(srt,1,ndim-nocc+1)
         norm = norm + (1d0-ss4)*0.10d0 
         fitdat(fitcount) = norm
         totmatch = totmatch + norm
         if(pr)write(*,'(2i4,6f11.4)') i,match,epsref(match),e(i),norm,-srt(1),-srt(2)
      enddo 

      if(pr)then
         write(*,*) 'total MO match score incl virt:', totmatch  
         write(*,'('' gap (eV)  DFT vs. gTB        : '',2f9.5)') (epsref(nocc+1)-epsref(nocc))*au2ev,(e(nocc+1) - e(nocc))*au2ev
      endif

      end
