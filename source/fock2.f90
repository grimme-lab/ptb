subroutine fock2(n,ndim,at,xyz,z,rab,cn,S,SS,Vecp,Hdiag,&
                 scfpar,S1,S2,psh,pin,P,Hmat)
   use iso_fortran_env, only : wp => real64
   use bascom
   use parcom
   use com
   implicit none 
!! ------------------------------------------------------------------------
!  Input
!! ------------------------------------------------------------------------
   integer, intent(in)    :: n                     ! number of atoms 
   integer, intent(in)    :: ndim                  ! number of AOs       
   integer, intent(in)    :: at(n)                 ! ordinal number of atoms
   real(wp),intent(in)    :: xyz(3,n)              ! coordinates
   real(wp),intent(in)    :: z(n)                  ! nuclear charges          
   real(wp),intent(in)    :: rab(n*(n+1)/2)        ! distances  
   real(wp),intent(in)    :: cn(n)                 ! CN           
   real(wp),intent(in)    :: S(ndim*(ndim+1)/2)    ! exact overlap maxtrix in SAO
   real(wp),intent(in)    :: SS(ndim*(ndim+1)/2)   ! scaled overlap maxtrix in SAO
   real(wp),intent(in)    :: P   (ndim*(ndim+1)/2) ! density matrix from iter. 2
   real(wp),intent(in)    :: Vecp(ndim*(ndim+1)/2) ! ECP ints
   real(wp),intent(in)    :: Hdiag(ndim)           ! diagonal of H0
   real(wp),intent(in)    :: scfpar(8)             ! parameters
   real(wp),intent(in)    :: psh(10,n)             ! shell populations from iter. 2
   real(wp),intent(in)    :: pin(n)                ! atom      "         "   "    "
   real*4  ,intent(in)    :: S1(ndim,ndim)         ! ML trafo
   real*4  ,intent(in)    :: S2(ndim,ndim)         ! "   "          

!! ------------------------------------------------------------------------
!  Output
!! ------------------------------------------------------------------------
   real(wp),intent(inout) :: Hmat(ndim*(ndim+1)/2) ! Hamiltonian matrix

!  local
   logical  :: fail
   integer  :: i,j,k,l,ish,ati,atj,ia,ib,jsh,ii,jj,lin,ij,li,lj,iter,iish,jjsh,mode
!  real(wp),parameter :: cok   = 0.95_wp  ! OK mixing in 1. iter, not as important as in second iter
!  real(wp),parameter :: cmn   = 1_wp-cok     
   real(wp),parameter :: au2ev = 27.2113957_wp
   real(wp) :: r,tmp,pol,hi,hj,hij,xk,t8,t9,qa,qb,keav,eh1,dmp,tmp2,ssh
   real(wp) :: gq(n), geff(n), xab(n*(n+1)/2), pa(n)
   real(wp) :: ves(nsh), gab(nsh,nsh)  

   pa = z - pin  ! output are populations but here q is used for convenience
   call setgab  (n,at,rab,pa,scfpar(5),gab)  ! the gab contain q as higher order effect on Ves
   call setespot(n,at,psh,gab,ves) 
   ves = ves * 0.5_wp

! H0 +  third-order (atomic charge exists in 1. AND 2. iter)
   do i=1, n
      geff(i) = pa(i)**2*shell_cnf4(1,at(i)) ! geff is temp.
   enddo

   Hmat = 0_wp
   ij = 0
   do i=1,ndim
      ia = aoat(i)
      ati= at(ia)
      ish= shell2ao(i)
      li = bas_lsh(ish,ati)
      hi = Hdiag(i)
      do j=1,i  
         ij = ij + 1
         ib = aoat(j)
         r  = rab(lin(ia,ib))
         if(r.gt.50d0) cycle
         hj = Hdiag(j)
         hij= hi+hj
         ssh= hij * SS(ij)
         atj= at(ib)
         jsh= shell2ao(j)
         if(ia.ne.ib) then            ! different atoms
            lj  = bas_lsh(jsh,atj)
            xk  = (shell_cnf4(2,ati)+shell_cnf4(2,atj)) 
            pol = ((hi-hj)/hij)**2
            keav= 0.5_wp*(shell_cnf2(9,ati)+dble(li)*scfpar(2) + shell_cnf2(9,atj)+dble(lj)*scfpar(2))
            tmp = ssh * keav * (1_wp-pol*scfpar(1)) * (1_wp+xk/r) ! fit yields same values for iter1,2
         else                         ! same atoms
            if(ish.ne.jsh) then       ! s-s', p-p', d-d' off-diagonal, li=lj because S=0 otherwise
               tmp2= shell_cnf4(4+li,ati) 
               tmp = ssh * tmp2 + shell_cnf3(9,ati)* tmp2 * hij * SS(ij)**2 ! second term only for more than 2 shells of same l
            else
               tmp = ssh
            endif
         endif
!                                               third order diagonal
         Hmat(ij) = Hmat(ij) + tmp + Vecp(ij) - S(ij)*(geff(ia)+geff(ib))
      enddo
   enddo

! H1
!   for +U LR damping
    k = 0
    do i=1,n
      gq(i) = 1_wp-(shell_xi(9,at(i))*pa(i)+shell_xi(10,at(i))*pa(i)**2) ! gq is temp., important charge scaling
      hi = shell_cnf3(10,at(i)) + (cn(i)-avcn(at(i)))*shell_resp(10,at(i),1)
      do j=1,i
         k = k + 1
         r = hi + shell_cnf3(10,at(j)) + (cn(j)-avcn(at(j)))*shell_resp(10,at(j),1)
         t8= (rab(k)-r)/r
         xab(k) = 0.5_wp*(1_wp+erf(-1.8_wp*t8)) ! paramter not very important
      enddo
    enddo
    call calcpauli2(n,ndim,at,psh,S,Hdiag,Hmat) ! add valence X correction based on three-index ECP formula, now using psh
    k = 0
    do i=1,ndim
      ia = aoat(i)
      ati= at(ia)
      ish= shell2ao(i)
      iish=shmap(ish,ia)
      hi = shell_cnf2(ish,ati)*gq(ia)
      do j=1,i-1
         k  = k + 1
         ib = aoat(j)
         atj= at(ib)
         jsh= shell2ao(j)
         jjsh=shmap(jsh,ib)
         hj = shell_cnf2(jsh,atj)*gq(ib)
!                            this part is INDO two-c like         shell ES
         Hmat(k) = Hmat(k) + P(k) * (hi + hj) * xab(lin(ib,ia)) - S(k)*(ves(iish)+ves(jjsh)) 
      enddo
      k = k + 1
      Hmat(k) = Hmat(k) + 2d0*P(k) * shell_xi(8,ati) * hi * xab(lin(ia,ia)) - 2d0*S(k)*ves(iish) ! scaled diag
    enddo

end

