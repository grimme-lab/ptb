subroutine onescf(n,ndim,nel,nopen,homo,at,shmap,rab,cn,S,SS,Hmat,Hdiag,focc,&
                  eT,scfpar,ves0,psh,pa,P)                      
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
   integer, intent(in)    :: nel                   ! number of electrons 
   integer, intent(in)    :: nopen                 ! number of open shells
   integer, intent(in)    :: homo                  ! as the name says...
   integer, intent(in)    :: at(n)                 ! ordinal number of atoms
   integer, intent(in)    :: shmap(10,n)           ! 
   real(wp),intent(in)    :: rab(n*(n+1)/2)        ! distances  
   real(wp),intent(in)    :: cn(n)                 ! CN           
   real(wp),intent(in)    :: S(ndim*(ndim+1)/2)    ! exact overlap maxtrix in SAO
   real(wp),intent(in)    :: SS(ndim*(ndim+1)/2)   ! scaled overlap maxtrix in SAO
   real(wp),intent(inout) :: Hmat(ndim*(ndim+1)/2) ! Vecp + field initialized
   real(wp),intent(in)    :: Hdiag(ndim)           ! diagonal of H0
   real(wp),intent(in)    :: focc (ndim)           ! fractional occ.
   real(wp),intent(in)    :: eT                    ! el. temp.
   real(wp),intent(in)    :: scfpar(8)             ! parameters
   real(wp),intent(in)    :: ves0(nsh)             ! ES potential field free      
   real(wp),intent(in)    :: psh(10,n)             ! shell populations with field
   real(wp),intent(in)    :: pa(n)                 ! atom      "         "   "

!! ------------------------------------------------------------------------
!  Output
!! ------------------------------------------------------------------------
   real(wp),intent(inout)   :: P (ndim*(ndim+1)/2) ! density matrix

!  local
   logical  :: fail
   integer  :: i,j,k,l,ish,ati,atj,ia,ib,jsh,ii,jj,lin,ij,li,iish,jjsh,mode
   real(wp) :: r,tmp,pol,hi,hj,hij,xk,t8,t9,qa,qb,keav,eh1,dmp,tmp2,ssh,eel
   real(wp) :: vi,vj         
   real(wp) :: gq(n),geff(n),ves(nsh),eps(ndim),xab(nsh,nsh)

   call setgab  (n,at,rab,pa,scfpar(5),xab)  ! the gab contain q as higher order effect on Ves
   call setespot(n,at,psh,xab,ves) 
   ves = ves * 0.5_wp

! H0 +  third-order (atomic charge exists in 1. AND 2. iter)
   do i=1, n
      geff(i) = pa(i)**2*shell_cnf4(1,at(i)) ! geff is temp.
   enddo

   ij = 0
   do i=1,ndim
      ia = aoat(i)
      ati= at(ia)
      ish= shell2ao(i)
      li = bas_lsh(ish,ati)+1
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
         if(ia.ne.ib) then            ! different atoms
            xk  = (shell_cnf4(2,ati)+shell_cnf4(2,atj)) 
            pol = ((hi-hj)/hij)**2
            keav= 0.5_wp*(shell_cnf2(9,ati) + shell_cnf2(9,atj))
            tmp = ssh * keav * (1_wp-pol*scfpar(1)) * (1_wp+xk/r) ! fit yields same values for iter1,2
         else                         ! same atoms
            jsh = shell2ao(j)
            if(ish.ne.jsh) then       ! s-s', p-p', d-d' off-diagonal, li=lj because S=0 otherwise
               tmp2= shell_cnf4(3+li,ati) 
               tmp = ssh * tmp2 + shell_cnf3(9,ati)* tmp2 * hij * SS(ij)**2
            else
               tmp = ssh 
            endif
         endif
!                                               third order diagonal
         Hmat(ij) = Hmat(ij) + tmp - S(ij)*(geff(ia)+geff(ib))
      enddo
   enddo

! H1
    k = 0
    do i=1,n
      gq(i) = 1_wp-(shell_xi(9,at(i))*pa(i)+shell_xi(10,at(i))*pa(i)**2) 
      hi = shell_cnf3(10,at(i)) + (cn(i)-avcn(at(i)))*scfpar(7)
      do j=1,i
         k = k + 1
         r = hi + shell_cnf3(10,at(j)) + (cn(j)-avcn(at(j)))*scfpar(7)
         t8= (rab(k)-r)/r
         xab(j,i) = 0.5_wp*(1_wp+erf(-2_wp*t8)) 
         xab(i,j) = xab(j,i)
      enddo
    enddo
    call calcpauli2(n,ndim,at,psh,S,Hdiag,Hmat) 

    k = 0
    do i=1,ndim
      ia = aoat(i)
      ati= at(ia)
      ish= shell2ao(i)
      iish=shmap(ish,ia)
      hi = shell_cnf2(ish,ati)*gq(ia)*expscal(3,10,ati)                  ! adapted +U scaling
      vi = ves(iish)*expscal(3,9,ati)+ves0(iish)*(1_wp-expscal(3,9,ati)) ! mixing of field perturbed and field free ES potential 
      do j=1,i                                                           ! with element wise parameter
         k  = k + 1
         ib = aoat(j)
         atj= at(ib)
         jsh= shell2ao(j)
         jjsh=shmap(jsh,ib)
         vj = ves(jjsh)*expscal(3,9,atj)+ves0(jjsh)*(1_wp-expscal(3,9,atj))
         hj = shell_cnf2(jsh,atj)*gq(ib)*expscal(3,10,atj)
!                            this part is INDO two-c like         shell ES
         Hmat(k) = Hmat(k) + P(k) * (hi + hj) * xab(ib,ia) - S(k)*(vi+vj) 
      enddo
    enddo

   call solve2 (mode,ndim,nel,nopen,homo,eT,focc,Hmat,S,P,eps,eel,fail) 
   if(fail) stop 'diag error onescf'

end

