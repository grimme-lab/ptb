module sm
   use metrics
   use gtb_accuracy, only: ik, sp, dp, i4, wp
   use gtb_lapack_eig, only : la_syevd,la_sygvd
   use gtb_la, only: la_gemm, la_symm
   use timer, only : tTimer
   use omp_lib, only: omp_get_wtime
   implicit none
  type :: eigdecomp_type
    integer :: ndim
    real(wp),allocatable :: U(:,:)
    real(wp),allocatable :: e(:)
    real(wp),allocatable :: fe(:)
    real(wp),allocatable :: occ(:)
    integer,allocatable :: map(:)
    integer,allocatable :: imap(:)
  endtype

contains
         
   subroutine initseed()
        INTEGER :: i, n
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed
        real(8) :: svar
        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))

        seed = 42 * (/ (i - 1, i = 1, n) /)
        CALL RANDOM_SEED(PUT = seed)
        DEALLOCATE(seed)
        call random_number(svar)
        end subroutine

   subroutine sm_stupid_simple(ndim,nelref,Hvec,Svec,Pvec,n,xyz,aoat,submatrix_columns,submatrix_mode)
     integer,intent(in)     :: nelref
     integer,intent(in)     :: ndim
     real(wp),intent(in)    :: Hvec(:)
     real(wp),intent(in)    :: Svec(:)
     real(wp),intent(inout) :: Pvec(:)
     integer,intent(in)     :: n
     real(wp),intent(in)    :: xyz(3,n)
     integer,intent(in)     :: aoat(ndim)
     integer,intent(in)     :: submatrix_columns
     integer,intent(in)     :: submatrix_mode

     real(wp),allocatable :: H(:,:)
     real(wp),allocatable :: S(:,:)
     real(wp),allocatable :: P(:,:)
     real(wp),allocatable :: Hsm(:,:)
     real(wp),allocatable :: Ssm(:,:)
     real(wp),allocatable :: sqrtinvSsm(:,:)
     real(wp),allocatable :: Psm(:,:)
     real(wp) :: effort_full, effort_sms, effort_min
     integer,allocatable :: comb(:,:)
     integer,allocatable :: ncomb(:)
     integer,allocatable :: itmp(:)
     integer,allocatable :: tmp(:)
     integer,allocatable :: st(:)
     
     integer,allocatable :: kmeans_z(:)
     integer,allocatable :: kmeans_z_opt(:)
     real(wp),allocatable :: kmeans_c(:,:)
     real(wp),allocatable :: kmeans_work(:)
     real(wp),allocatable :: efforts(:)
     type(eigdecomp_type),allocatable :: eigdecomp(:)

     integer :: ism,i,j,k,ind,smdim,muiter,nsm,sumsmdims,nsm_opt
     real(wp) :: mua,mub,mu,nel,t0,t1
     type(tTimer) :: timer_sm
     integer(ik) :: info
     logical :: finaliter
     logical :: firstiter

     print*,"Stupid simple implementation of the submatrix method:",submatrix_columns,submatrix_mode
     call timer_sm%new(20)
     call timer_sm%click(1, 'alloc')

     allocate(H(ndim,ndim))
     allocate(S(ndim,ndim))
     allocate(P(ndim,ndim))
     allocate(itmp(ndim))
     allocate(tmp(ndim))
     allocate(st(ndim))
     
     call blowsym(ndim,Hvec,H)
     call blowsym(ndim,Svec,S)

     call timer_sm%click(1)
     call timer_sm%click(2, 'submatrix combination')

     nsm=0
     if(submatrix_columns.eq.1)then
       !columns as submatrices
       nsm=ndim
       allocate(comb(nsm,1))
       allocate(ncomb(nsm))
       ncomb(:)=1
       do i=1,ndim
         comb(i,1)=i
       enddo
     elseif(submatrix_columns.eq.2)then
       !atoms as submatrices
       allocate(comb(ndim,ndim))
       allocate(ncomb(ndim))
       ism=1
       ncomb(ism)=1
       comb(ism,ncomb(ism))=1

       do i=2,ndim
         if(aoat(i).eq.aoat(i-1))then
            ncomb(ism)=ncomb(ism)+1
            comb(ism,ncomb(ism))=i
         else
           ism=ism+1
           ncomb(ism)=1
           comb(ism,ncomb(ism))=i
         endif
       enddo
       nsm=ism
     elseif(submatrix_columns.eq.3)then
       !full matrix
       nsm=1
       allocate(comb(nsm,ndim))
       allocate(ncomb(nsm))
       ncomb(1)=ndim
       do i=1,ndim
         comb(1,i)=i
       enddo
     elseif(submatrix_columns.eq.4)then
       !coordinate-based clustering with simple kmeans (alternative: metis)
       effort_min=huge(effort_min)
       allocate(kmeans_z(n))
       allocate(kmeans_z_opt(n))
       allocate(kmeans_work(n))
       allocate(kmeans_C(3,n))
       allocate(comb(n,ndim))
       allocate(ncomb(n))
       !$OMP PARALLEL DO private(t0,t1,kmeans_c,kmeans_z,kmeans_work,effort_sms,ncomb,comb)
       do nsm=1,min(64,n)
         t0=omp_get_wtime()
         call clustering_kmeans(n,ndim,xyz,kmeans_c,kmeans_z,kmeans_work,H,S,aoat,nsm,ncomb,comb,effort_sms)
         t1=omp_get_wtime()
         print*,"kmeans",nsm,effort_sms,t1-t0
         if(effort_sms.lt.effort_min)then
           !$OMP Critical
           effort_min=effort_sms
           kmeans_z_opt(:)=kmeans_z(:)
           nsm_opt=nsm
           !$OMP end Critical
         endif
       enddo
       !$OMP end PARALLEL DO 
       nsm=nsm_opt
       call clustering_kmeans(n,ndim,xyz,kmeans_c,kmeans_z,kmeans_work,H,S,aoat,nsm,ncomb,comb,effort_sms)

       deallocate(kmeans_z)
       deallocate(kmeans_work)
       deallocate(kmeans_C)
     else
       print*,"not implemented"
       stop
     endif
     call timer_sm%click(2)
     call timer_sm%click(3, 'submatrix stats')
     !collect submatrix infos
     effort_full=real(ndim,kind=wp)**2.5D0
     effort_sms=0
     sumsmdims=0
     do ism=1,nsm
      ind=0
      do i=1,ndim
        do j=1,ncomb(ism)
         if(H(i,comb(ism,j)).gt.0.or.S(i,comb(ism,j)).gt.0.or.i.eq.comb(ism,j))then
           ind=ind+1
           exit
         endif
        enddo
      enddo
      smdim=ind
      sumsmdims=sumsmdims+smdim
      effort_sms=effort_sms+real(smdim,kind=wp)**2.5D0
     enddo
     print*,"number of submatrices=",nsm
     print*,"average sm dim=",sumsmdims/nsm
     print*,"effort for full matrix ops=",effort_full
     print*,"sums of efforts for submatrices=",effort_sms
     print*,"ratio",effort_sms/effort_full
     call timer_sm%click(3)
     
     call timer_sm%click(4, 'mu iter')
     mua=-10
     mub=10
     firstiter=.true.
     finaliter=.false.
     allocate(eigdecomp(nsm))

     do muiter=1,100
       mu=(mua+mub)/2

       nel=0
       if(finaliter)then
         P(:,:)=0
       endif

       !submatrix method
       do ism=1,nsm
        if(firstiter)then
          call timer_sm%click(5, 'submatrix construction')
          !collect submatrix infos
          st(:)=0
          !$OMP PARALLEL DO Private(j)
          do i=1,ndim
            do j=1,ncomb(ism)
             if(H(i,comb(ism,j)).gt.0.or.S(i,comb(ism,j)).gt.0.or.i.eq.comb(ism,j))then
               st(i)=1
               exit
             endif
            enddo
          enddo
          !$OMP end PARALLEL DO 
          smdim=0
          do i=1,ndim
            if(st(i).eq.1)then
               tmp(smdim+1)=i 
               itmp(i)=smdim+1
               smdim=smdim+1
            endif
          enddo

          allocate(eigdecomp(ism)%map(smdim))
          allocate(eigdecomp(ism)%imap(ndim))
          eigdecomp(ism)%ndim=smdim
          eigdecomp(ism)%map(1:smdim)=tmp(1:smdim)
          eigdecomp(ism)%imap(1:ndim)=itmp(1:ndim)

          !build submatrices
          allocate(Hsm(smdim,smdim))
          allocate(Ssm(smdim,smdim))
          allocate(Psm(smdim,smdim))
          !$OMP PARALLEL DO private(j)
          do i=1,smdim
            do j=1,smdim
              Hsm(i,j)=H(eigdecomp(ism)%map(i),eigdecomp(ism)%map(j))    
              Ssm(i,j)=S(eigdecomp(ism)%map(i),eigdecomp(ism)%map(j))    
            enddo
          enddo
          call timer_sm%click(5)
          call timer_sm%click(6, 'submatrix matrix function')
          if(submatrix_mode.eq.1)then
            print*,"submatrix_mode.eq.1 not implemented"
            stop
!            allocate(sqrtinvSsm(smdim,smdim))
!            call timer_sm%click(4, 'submatrix orthogonalization')
!            !Ortho
!            call matrix_root(smdim,Ssm,sqrtinvSsm)
!    !        Hsm(:,:)=matmul(matmul(sqrtinvSsm,Hsm),sqrtinvSsm)
!            call la_gemm(sqrtinvSsm, Hsm, Psm)
!            call la_gemm(Psm,sqrtinvSsm, Hsm)
!            
!            call timer_sm%click(4)
!            call timer_sm%click(5, 'submatrix sign')
!            
!            do i=1,smdim
!              Hsm(i,i)=Hsm(i,i)-mu
!            enddo
!            !sign function
!            call matrix_sign(smdim,Hsm,Ssm)
!
!            call timer_sm%click(5)
!            call timer_sm%click(6, 'submatrix nel')
!
!            !electron number
!            do j=1,ncomb(ism)
!              k=imap(comb(ism,j))
!              nel=nel+(1.0D0-Ssm(k,k))
!            enddo
!
!            call timer_sm%click(6)
!            call timer_sm%click(7, 'submatrix P')
!            !compute Psm
!            Ssm(:,:)=-Ssm(:,:)
!            do i=1,smdim
!              Ssm(i,i)=1.0D0+Ssm(i,i)
!            enddo
!    !        Psm(:,:)=matmul(sqrtinvSsm,matmul(Ssm,sqrtinvSsm))
!            call la_gemm(sqrtinvSsm, Ssm, Hsm)
!            call la_gemm(Hsm,sqrtinvSsm, Psm)
!            deallocate(sqrtinvSsm)
!            call timer_sm%click(7)
          elseif(submatrix_mode.eq.2)then
            call timer_sm%click(7, 'submatrix eigendecomp')
            allocate(eigdecomp(ism)%U(smdim,smdim))
            allocate(eigdecomp(ism)%e(smdim))
            allocate(eigdecomp(ism)%fe(smdim))
            allocate(eigdecomp(ism)%occ(smdim))
            eigdecomp(ism)%U(:,:)=Hsm(:,:)
            Psm(:,:)=Ssm(:,:)
            t0=omp_get_wtime()
            call la_sygvd(eigdecomp(ism)%U,Psm,eigdecomp(ism)%e,INFO)
            t1=omp_get_wtime()
            print*,"sm",ism,smdim,t1-t0
            call timer_sm%click(7)

            call timer_sm%click(8, 'submatrix fe')
            call la_gemm(Ssm,eigdecomp(ism)%U, Psm)
            eigdecomp(ism)%fe(:)=0
            do k=1,ncomb(ism)
              i=eigdecomp(ism)%imap(comb(ism,k))
              eigdecomp(ism)%fe(:)=eigdecomp(ism)%fe(:)+Psm(i,:)*eigdecomp(ism)%U(i,:)
            enddo
            call timer_sm%click(8)
          else
            print*,"not implemented submatrix_mode",submatrix_mode
            stop
          endif
          call timer_sm%click(6)
          deallocate(Hsm)
          deallocate(Ssm)
          deallocate(Psm)
        endif 
        call timer_sm%click(9, 'submatrix occ')
        do i=1,eigdecomp(ism)%ndim
          if(eigdecomp(ism)%e(i).lt.mu)then
            eigdecomp(ism)%occ(i)=2
          elseif(eigdecomp(ism)%e(i).eq.mu)then
            eigdecomp(ism)%occ(i)=1
          else
            eigdecomp(ism)%occ(i)=0
          endif
        enddo
        call timer_sm%click(9)
        
        call timer_sm%click(10, 'submatrix Nel')
        !call dmat(smdim,eigdecomp(ism)%occ,eigdecomp(ism)%U,Psm)
        !call la_gemm(Psm, Ssm, Hsm)
        !do j=1,ncomb(ism)
        !  k=imap(comb(ism,j))
        !  nel=nel+Hsm(k,k)
        !enddo
        nel=nel+sum(eigdecomp(ism)%fe(:)*eigdecomp(ism)%occ(:))
        call timer_sm%click(10)
        
        if(finaliter)then
          call timer_sm%click(11, 'submatrix final iter')
          smdim=eigdecomp(ism)%ndim
          allocate(Psm(smdim,smdim))
          do i=1,smdim
            do j=1,smdim
              Psm(i,j)=P(eigdecomp(ism)%map(i),eigdecomp(ism)%map(j))    
            enddo
          enddo
          !put results in P
          call dmat(smdim,eigdecomp(ism)%occ,eigdecomp(ism)%U,Psm)
          do j=1,ncomb(ism)
            k=eigdecomp(ism)%imap(comb(ism,j))
            do i=1,smdim
              P(eigdecomp(ism)%map(i),comb(ism,j))=Psm(i,k)
            enddo
          enddo
          deallocate(Psm)
          call timer_sm%click(11)
        endif
       
       enddo !submatrix iteration
       if(finaliter)then
        exit
       endif
       
       print*,"muiter",muiter,"mu",mu,nel,"of",nelref
       if(abs(mua-mub).lt.1e-3.or.abs(nel-nelref).lt.1e-5)then
         print*,"iteration of chemical potential converged"
         finaliter=.true.
         print*,"final iteration to get density matrix"
       endif
       if(.not.finaliter)then
         if(nel>nelref)then
           mub=mu
         else
           mua=mu
         endif
       endif
       firstiter=.false.
     enddo !mu iteration
     call timer_sm%click(4)

     call timer_sm%click(12, 'submatrix dealloc')
     call packsym(ndim,P,Pvec)
     deallocate(H)
     deallocate(S)
     deallocate(P)
     deallocate(ncomb)
     deallocate(comb)
     if(allocated(eigdecomp))then
       do ism=1,nsm
         if(allocated(eigdecomp(ism)%U))then
           deallocate(eigdecomp(ism)%U)
           deallocate(eigdecomp(ism)%e)
           deallocate(eigdecomp(ism)%fe)
           deallocate(eigdecomp(ism)%occ)
           deallocate(eigdecomp(ism)%map)
           deallocate(eigdecomp(ism)%imap)
         endif
       enddo
       deallocate(eigdecomp)
       deallocate(st,itmp,tmp)
     endif
     call timer_sm%click(12)
     call timer_sm%finalize('Submatrix')
   end subroutine
   subroutine matrix_sign(ndim,A,S)
     integer,intent(in) :: ndim
     real(wp),intent(inout) :: A(ndim,ndim)
     real(wp),intent(inout) :: S(ndim,ndim)
     real(wp) :: e(ndim)
     real(wp) :: D(ndim,ndim)
     integer(ik) :: info
     integer :: i

     e(:)=1
     call la_syevd(A,e,INFO)
     D(:,:)=0
     do i=1,ndim
       if(e(i).gt.0)then
         D(i,i)=1
       elseif(e(i).lt.0)then
         D(i,i)=-1
       else
         D(i,i)=0
       endif
     enddo
     call la_gemm(A, D, S)
     D(:,:)=S(:,:)
     call la_gemm(D, A, S,transb='T')
!     S=matmul(A,matmul(D,transpose(A)))
     
   end subroutine
   subroutine matrix_root(ndim,A,sqrtAinv)
     integer,intent(in) :: ndim
     real(wp),intent(inout) :: A(ndim,ndim)
     real(wp),intent(inout) :: sqrtAinv(ndim,ndim)
     real(wp) :: e(ndim)
     real(wp) :: D(ndim,ndim)
     integer(ik) :: info
     integer :: i

     e(:)=1
     call la_syevd(A,e,INFO)
     !print*,"matrix_root: eigenvalues",info,e
     D(:,:)=0
     do i=1,ndim
       D(i,i)=1.0D0/sqrt(e(i))
     enddo
     call la_gemm(A, D, sqrtAinv)
     D(:,:)=sqrtAinv(:,:)
     call la_gemm(D, A, sqrtAinv,transb='T')
     !sqrtAinv=matmul(A,matmul(D,transpose(A)))
   end subroutine
   
   subroutine clustering_kmeans(n,ndim,xyz,kmeans_c,kmeans_z,kmeans_work,H,S,aoat,nsm,ncomb,comb,effort_sms)
     integer,intent(in) :: n,ndim
     real(wp),intent(in) :: xyz(3,n)
     real(wp),intent(inout) :: kmeans_c(3,n)
     integer,intent(inout) :: kmeans_z(n)
     real(wp),intent(inout) :: kmeans_work(n)
     real(wp),intent(in) :: H(ndim,ndim)
     real(wp),intent(in) :: S(ndim,ndim)
     integer,intent(in) :: aoat(ndim)
     integer,intent(in) :: nsm
     real(wp),intent(out) :: effort_sms
     integer,intent(inout) :: ncomb(:)
     integer,intent(inout) :: comb(:,:)
     integer :: i,ism,ind,j,smdim

     kmeans_C(:,:)=0
     kmeans_Z(:)=0
     call initseed()
     CALL KMPP(xyz, 3, n, nsm, kmeans_c, kmeans_Z, kmeans_WORK, I)
     if(i.ne.0)then
       print*,"kmeans failed",I
       stop
     endif
     !build combination information from clustering
     ncomb(:)=0
     do ism=1,nsm
       do i=1,n !atoms
         if(kmeans_z(i).eq.ism)then
           !atom is in submatrix
           do j=1,ndim !columns
             if(aoat(j).eq.i)then
               ncomb(ism)=ncomb(ism)+1
               comb(ism,ncomb(ism))=j
             endif
           enddo
         endif
       enddo
     enddo

     !evaluate effort
     effort_sms=0
     do ism=1,nsm
       ind=0
       do i=1,ndim
         do j=1,ncomb(ism)
           if(H(i,comb(ism,j)).gt.0.or.S(i,comb(ism,j)).gt.0.or.i.eq.comb(ism,j))then
              !map(ind+1)=i 
              ind=ind+1
              exit
           endif
         enddo
       enddo
       smdim=ind
       effort_sms=effort_sms+real(smdim,kind=wp)**2.5D0
     enddo
   end subroutine
end module
