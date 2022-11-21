
!------------------------------------------------------------------------------------------------------
! distributed atomic multipole moment interactions: all interactions up to r**-3
! energy evaluation
! nat         : # of atoms
! xyz(3,nat)  : cartesian coordinates
! qm(nat)     : atomic partial charges
! dipm(3,nat) : cumulative atomic dipole moments (x,y,z)
! qp(6,nat)   : traceless(!) cumulative atomic quadrupole moments (xx,xy,yy,xz,yz,zz)
! e           : E_AES
! local
! gab3,gab5   : damped R**-3 and R**-5 Coulomb laws, dimension: nat*(nat+1)/2
!               multiplication with numerator then leads to R**-2 and R**-3 decay, respectively
!------------------------------------------------------------------------------------------------------

      subroutine aniso_es(nat,at,rab,xyz,e)
      use iso_fortran_env, only : wp => real64
      use com
      use parcom
      use aescom, only: qm,dipm,qp
      implicit none          
      integer nat,at(nat)
      real*8 rab(nat*(nat+1)/2)
      real*8 xyz(3,nat)
      real*8 e           

      real*8 qp1(6),rr(3),dp1(3),rij(3)
      real*8 edd,e01,e02,e11,r0,r5,r3,r2,tt,tt3,q1,qs2,xk,damp
      real*8 ed,eq
      real*8 gab3(nat*(nat+1)/2),gab5(nat*(nat+1)/2)
      integer i,j,k,l,m,ki,kj,kl,lin

   real(wp),parameter :: rcov_mod(118) = 1.889725949_wp * [ & 
   & 0.17_wp,0.46_wp, & ! H changed
   & 1.20_wp,0.94_wp,0.77_wp,0.85_wp,0.68_wp,0.62_wp,0.56_wp,0.50_wp, & ! Li-Ne CNOF hand adjusted
   & 1.40_wp,1.25_wp,1.13_wp,1.14_wp,1.08_wp,1.25_wp,0.80_wp,0.96_wp, & ! Na-Ar, Si-Cl changed
   & 1.76_wp,1.54_wp, & ! K,Ca
   &                 1.33_wp,1.22_wp,1.21_wp,1.10_wp,1.07_wp, & ! Sc-
   &                 1.04_wp,1.00_wp,0.99_wp,1.01_wp,1.09_wp, & ! -Zn
   &                 1.12_wp,1.09_wp,1.15_wp,1.10_wp,1.14_wp,1.17_wp, & ! Ga-Kr
   & 1.89_wp,1.67_wp, & ! Rb,Sr
   &                 1.47_wp,1.39_wp,1.32_wp,1.24_wp,1.15_wp, & ! Y-
   &                 1.13_wp,1.13_wp,1.08_wp,1.15_wp,1.23_wp, & ! -Cd
   &                 1.28_wp,1.26_wp,1.26_wp,1.23_wp,1.32_wp,1.31_wp, & ! In-Xe
   & 2.09_wp,1.76_wp, & ! Cs,Ba
   &         1.62_wp,1.47_wp,1.58_wp,1.57_wp,1.56_wp,1.55_wp,1.51_wp, & ! La-Eu
   &         1.52_wp,1.51_wp,1.50_wp,1.49_wp,1.49_wp,1.48_wp,1.53_wp, & ! Gd-Yb
   &                 1.46_wp,1.37_wp,1.31_wp,1.23_wp,1.18_wp, & ! Lu-
   &                 1.16_wp,1.11_wp,1.12_wp,1.13_wp,1.32_wp, & ! -Hg
   &                 1.30_wp,1.30_wp,1.36_wp,1.31_wp,1.38_wp,1.42_wp, & ! Tl-Rn
   & 2.01_wp,1.81_wp, & ! Fr,Ra
   &         1.67_wp,1.58_wp,1.52_wp,1.53_wp,1.54_wp,1.55_wp,1.49_wp, & ! Ac-Am
   &         1.49_wp,1.51_wp,1.51_wp,1.48_wp,1.50_wp,1.56_wp,1.58_wp, & ! Cm-No
   &                 1.45_wp,1.41_wp,1.34_wp,1.29_wp,1.27_wp, & ! Lr-
   &                 1.21_wp,1.16_wp,1.15_wp,1.09_wp,1.22_wp, & ! -Cn
   &                 1.36_wp,1.43_wp,1.46_wp,1.58_wp,1.48_wp,1.57_wp ] ! Nh-Og

      k=0
      do i=1, nat
         do j=1,i-1
            k = k + 1
            r3= rab(k)**3
            r5= r3*rab(k)*rab(k)
!           xk= (ener_par5(8,at(i)) + ener_par5(8,at(j)))*0.5d0 ! small effect so take more or less standard radii
            r0=2.6d0*(rcov_mod(at(i))+rcov_mod(at(j)))          ! GLOB_PAR
            damp= 1d0 - 0.5d0*(erf(-1.5d0*(rab(k)-r0)/r0)+1d0)  ! GLOB_PAR
            gab3(k) = damp/r3
            gab5(k) = damp/r5
         enddo
         k=k+1
      enddo

      e01=0.0d0
      e02=0.0d0
      e11=0.0d0
      do i=2,nat 
         q1=qm(i)           ! charge scaling
         rr(1:3)=xyz(1:3,i)
         dp1(1:3)=dipm(1:3,i) 
         qp1(1:6)=qp(1:6,i)   
         do j=1,i-1      
            kj=lin(j,i)
            qs2=1d0 !                ! charge scaling
            rij(1:3)=xyz(1:3,j)-rr(1:3)
            r2=rab(kj)            
            ed=0.0d0
            eq=0.0d0
            edd=0.0d0
!           dipole - charge
            do k=1,3
               ed=ed+qm(j)*dp1(k)*rij(k)*qs2
               ed=ed-dipm(k,j)*q1*rij(k)
!              dip-dip & charge-qpole
               do l=1,3
                  kl=lin(l,k)
                  tt=rij(l)*rij(k)
                  tt3=3.0d0*tt
                  eq=eq+qm(j)*qp1(kl)*tt*qs2
                  eq=eq+qp(kl,j)*q1*tt
                  edd=edd-dipm(k,j)*dp1(l)*tt3 
               enddo
!              diagonal dip-dip term
               edd=edd+dipm(k,j)*dp1(k)*r2
            enddo
            e01=e01+ ed*gab3(kj)
            e02=e02+ eq*gab5(kj)
            e11=e11+edd*gab5(kj)
         enddo
      enddo

      e = e01 + e02 + e11

!     write(*,'(''d,q,dd'',3f9.5)')  e01,e02,e11

      end

!------------------------------------------------------------------------------------------------------
!      compute the cumulative atomic dipole and quadrupole moments via Mulliken population analysis
! input
! nat              : # of atoms
! nao              : # of spherical AOs (SAOs)
! xyz(3,nat)       : cartesian coordinates
! z(nat)           : nuclear charges      
! s                : exact overlap matrix
! p                : density matrix
! input on aescom 
!  pint            : dipole integral matrix, dimension 1-3,nao*(nao+1)/2
!                  : quadrupole integral matrix, dimension 4-9,nao*(nao+1)/2
! output on aescom 
! qm  (nat)        : Mulliken charges
! dipm(3,nat)      : cumulative atomic dipole moments (x,y,z)
! qp(6,nat)        : traceless(!) cumulative atomic quadrupole moments (xx,xy,yy,xz,yz,zz)
!------------------------------------------------------------------------------------------------------

subroutine camm(nat,nao,xyz,z,p,s)
      use bascom
      use aescom
      implicit none
      integer, intent(in) :: nao,nat
      real*8, intent(in) :: s(nao*(nao+1)/2),p(nao*(nao+1)/2)
      real*8, intent(in) :: xyz(3,nat)
      real*8, intent(in) :: z(nat)

      real*8 xk1,xl1,xk2,xl2,pij,tii,tjj
      real*8 pqm,pdmk,pdml,ps,ra(3)

      integer i,j,k,l,ii,jj,ij,kl,kj,lin

      qm  =0.0d0
      dipm=0.0d0
      qp  =0.0d0

      ij=0
      do i=1,nao
         ii=aoat(i)
         ra(1:3)=xyz(1:3,ii)
         do j=1,i-1
            ij=ij+1
            jj=aoat(j)
            pij=p(ij)
            ps=pij*s(ij)
            qm(ii)=qm(ii)+ps
            qm(jj)=qm(jj)+ps
            kl=0
            !  the qpint is stored as xx,yy,zz,xy,xz,yz (from integral routine)
            !  when doing the Mulliken population, we switch to lin-compatible sorting
            !  i,e. xx,xy,yy,xz,yz,zz
            do k=1,3
               xk1=ra(k)
               xk2=xyz(k,jj)   
               pdmk=pij*pint(k,ij)
               dipm(k,jj)=dipm(k,jj)+xk2*ps-pdmk
               dipm(k,ii)=dipm(k,ii)+xk1*ps-pdmk
               ! off-diagonal 
               do l=1,k-1
                  kl=kl+1
                  kj=k+l+1 
                  xl1=ra(l)
                  xl2=xyz(l,jj)
                  pdml=pij*pint(l,ij) 
                  pqm=pij*pint(3+kj,ij)
                  tii=pdmk*xl1+pdml*xk1-xl1*xk1*ps-pqm
                  tjj=pdmk*xl2+pdml*xk2-xl2*xk2*ps-pqm
                  qp(kl,jj)=qp(kl,jj)+tjj
                  qp(kl,ii)=qp(kl,ii)+tii
               enddo
               ! diagonal 
               kl=kl+1
               pqm=pij*pint(3+k,ij)
               tii=2.0d0*pdmk*xk1-xk1*xk1*ps-pqm
               tjj=2.0d0*pdmk*xk2-xk2*xk2*ps-pqm
               qp(kl,jj)=qp(kl,jj)+tjj
               qp(kl,ii)=qp(kl,ii)+tii
            enddo
         enddo
         ! diag in AO pair
         ij=ij+1
         pij=p(ij)
         ps=pij*s(ij)
         qm(ii)=qm(ii)+ps    
         kl=0
         do k=1,3
            xk1=ra(k)
            pdmk=pij*pint(k,ij)
            dipm(k,ii)=dipm(k,ii)+xk1*ps-pdmk
            ! off-diagonal 
            do l=1,k-1
               kl=kl+1
               kj=k+l+1 ! the qpint is stored as xx,yy,zz,xy,xz,yz (from integral routine)
               xl1=ra(l)
               pdml=pij*pint(l,ij)
               pqm=pij*pint(3+kj,ij)
               tii=pdmk*xl1+pdml*xk1-xl1*xk1*ps-pqm
               qp(kl,ii)=qp(kl,ii)+tii
            enddo
            !diagonal 
            kl=kl+1
            pqm=pij*pint(3+k,ij)
            tii=2.0d0*pdmk*xk1-xk1*xk1*ps-pqm
            qp(kl,ii)=qp(kl,ii)+tii
         enddo
      enddo
      ! remove trace
      do i=1,nat
         tii=qp(1,i)+qp(3,i)+qp(6,i)
         tii=0.50d0*tii
         qp(1:6,i)=1.50d0*qp(1:6,i)
         qp(1,i)=qp(1,i)-tii
         qp(3,i)=qp(3,i)-tii
         qp(6,i)=qp(6,i)-tii
         qm(i)  =z(i) - qm(i)  ! Mulliken charge
      enddo 

end subroutine camm   

!------------------------------------------------------------------------------------------------------
! add multipole ES potential to Fock matrix H
! data as energy routine
! NOT USED
!------------------------------------------------------------------------------------------------------


subroutine addaes(nat,ndim,vs,vd,vq,S,H)
   use iso_fortran_env, only : wp => real64
   use parcom, only: glob_par
   use bascom, only: aoat
   use aescom 
   implicit none
   integer, intent(in)   :: nat,ndim
   real(wp), intent(in)  :: vs(nat),vd(3,nat),vq(6,nat)
   real(wp), intent(in)  :: S(ndim*(ndim+1)/2)
   real(wp), intent(out) :: H(ndim*(ndim+1)/2)

   integer i,j,k,l,ii,jj
   real(wp) eh1,ff

   ff = 0.5d0 * glob_par(8) ! global scale factor, 1/2 taken out of sum term

   k = 0
   do i=1,ndim
      ii = aoat(i)
      do j=1,i
         jj = aoat(j)
         k = k + 1
         eh1=S(k)*(vs(ii)+vs(jj))                   ! charge
         H(k) = H(k) + eh1*ff   
         eh1 = 0
         do l=1,3
            eh1=eh1+pint(l,k)*(vd(l,ii)+vd(l,jj))   ! dipole
         enddo
         H(k) = H(k) + eh1*ff
         eh1 = 0
         do l=1,6
            eh1=eh1+pint(3+l,k)*(vq(l,ii)+vq(l,jj)) ! quadrupole
         enddo
         H(k) = H(k) + eh1*ff
      enddo
   enddo

end

!------------------------------------------------------------------------------------------------------
! precompute atomic multipole ES potential contributions vs,vd,vq to Fock matrix H
! data as energy routine
!------------------------------------------------------------------------------------------------------

subroutine setvsdq(nat,at,xyz,rab,vs,vd,vq)
   use iso_fortran_env, only : wp => real64
   use aescom 
   use parcom 
   use com 
   implicit none
   integer, intent(in)  :: nat,at(nat)
   real(wp), intent(in) :: xyz(3,nat)
   real(wp), intent(in) :: rab(nat*(nat+1)/2) 
   real(wp), intent(out) :: vs(nat),vd(3,nat),vq(6,nat)

   real(wp) ra(3),dra(3),rb(3),stmp,dum3a,dum5a,t1a,t2a,t3a,t4a,r2a
   real(wp) r2ab,t1b,t2b,t3b,t4b,dum3b,dum5b,dtmp(3),qtmp(6),g3,g5
   real(wp) qs1,qs2,damp,r3,r5,r0
   integer i,j,k,l1,l2,ll,m,mx,ki,kj
   integer lin

   vs = 0.0_wp
   vd = 0.0_wp
   vq = 0.0_wp
   ! set up overlap proportional potential
   do i = 1,nat
      ra(1:3) = xyz(1:3,i)
      stmp = 0.0_wp
      dtmp = 0.0_wp
      qtmp = 0.0_wp
      do j = 1,nat
         k  = lin(j,i)
         r3 = rab(k)**3
         r5 = r3*rab(k)*rab(k)
         r0 = glob_par(19)*(rcov(at(i))+rcov(at(j)))
         damp=1_wp - 0.5_wp*(erf(-1_wp*(rab(k)-r0)/r0)+1_wp)  
         g3 = damp/r3  
         g5 = damp/r5   
         rb(1:3) = xyz(1:3,j)
         dra(1:3) = ra(1:3)-rb(1:3)
         dum3a = 0.0_wp ! collect gab3 dependent terms
         dum5a = 0.0_wp ! collect gab5 dependent terms
         r2a = 0.0_wp
         r2ab = 0.0_wp
         t1a = 0.0_wp
         t2a = 0.0_wp
         t3a = 0.0_wp
         t4a = 0.0_wp
         ll = 0
         do l1 = 1,3
            ! potential from dipoles
            r2a = r2a+ra(l1)*ra(l1)      ! R_C * R_C
            r2ab = r2ab+dra(l1)*dra(l1)  ! R_AC * R_AC
            t1a = t1a+ra(l1)*dra(l1)     ! R_C * R_AC  : for dip-q (q-shift) and dip-dip (q-shift)
            t2a = t2a+dipm(l1,j)*dra(l1) ! mu_A * R_AC : for q-dip and dip-dip (q-shift)
            t3a = t3a+ra(l1)*dipm(l1,j)  ! R_C * mu_A  : for diag. dip-dip (q-shift)
            t4a = t4a+dra(l1)*dra(l1)*ra(l1)*ra(l1) ! (R_C o R_AC)**"2(square of Hadamard product) :
            ! results from trace remove from q-pole (q-shift)
            do l2 = 1,3
               ll = lin(l1,l2)
               ! potential from quadrupoles
               dum5a = dum5a-qp(ll,j)*dra(l1)*dra(l2) &
                  & -1.50_wp*qm(j)*dra(l1)*dra(l2)*ra(l1)*ra(l2)
               if(l2.ge.l1) cycle
               ki = l1+l2+1
               qtmp(ki) = qtmp(ki)-3.0_wp*qm(j)*g5*dra(l2)*dra(l1)
            enddo
            qtmp(l1) = qtmp(l1)-1.50_wp*qm(j)*g5*dra(l1)*dra(l1)
         enddo
         !
         ! set up S-dependent potential
         dum3a = -t1a*qm(j)-t2a ! dip-q (q-shift) and q-dip
         dum5a = dum5a+t3a*r2ab-3.0_wp*t1a*t2a & !dip-dip (q-shift terms)
            & +0.50_wp*qm(j)*r2a*r2ab !qpole-q (q-shift, trace removal)
         stmp = stmp+dum5a*g5+dum3a*g3
         do l1 = 1,3
            dum3a = dra(l1)*qm(j)
            dum5a = 3.0_wp*dra(l1)*t2a &           ! dipint-dip
               & -r2ab*dipm(l1,j) &            ! dipint-dip (diagonal)
               & -qm(j)*r2ab*ra(l1) &           ! qpole-q (dipint-shift, trace removal)
               & +3.0_wp*qm(j)*dra(l1)*t1a   ! qpole-q (dipint-shift)
            dtmp(l1) = dtmp(l1)+dum3a*g3+dum5a*g5
            qtmp(l1) = qtmp(l1)+0.50_wp*r2ab*qm(j)*g5 !remove trace term
         enddo
      enddo
      vs(i) = stmp                       ! q terms
      vd(1:3,i) = dtmp(1:3)              ! dipints from atom i
      vq(1:6,i) = qtmp(1:6)              ! qpints from atom i
      ! --- CT correction terms
      qs1 = 2.0_wp !aesData%dipKernel(at(i))*2.0_wp
      qs2 = 6.0_wp !aesData%quadKernel(at(i))*6.0_wp ! qpole pot prefactors
      t3a = 0.0_wp
      t2a = 0.0_wp
      do l1 = 1,3
         ! potential from dipoles
         t3a = t3a+ra(l1)*dipm(l1,i)*qs1  ! R_C * mu_C  : for diag. dip-dip
         vd(l1,i) = vd(l1,i)-qs1*dipm(l1,i)
         do l2 = 1,l1-1
            ! potential from quadrupoles
            ll = lin(l1,l2)
            ki = l1+l2+1
            vq(ki,i) = vq(ki,i)-qp(ll,i)*qs2
            t3a = t3a-ra(l1)*ra(l2)*qp(ll,i)*qs2
            vd(l1,i) = vd(l1,i)+ra(l2)*qp(ll,i)*qs2
            vd(l2,i) = vd(l2,i)+ra(l1)*qp(ll,i)*qs2
         enddo
         ! diagonal
         ll = lin(l1,l1)
         vq(l1,i) = vq(l1,i)-qp(ll,i)*qs2*0.50_wp
         t3a = t3a-ra(l1)*ra(l1)*qp(ll,i)*qs2*0.50_wp
         vd(l1,i) = vd(l1,i)+ra(l1)*qp(ll,i)*qs2
         ! collect trace removal terms
         t2a = t2a+qp(ll,i)
      enddo
      vs(i) = vs(i)+t3a
      ! trace removal
!     t2a = t2a*aesData%quadKernel(at(i))
      do l1 = 1,3
         vq(l1,i) = vq(l1,i)+t2a
         vd(l1,i) = vd(l1,i)-2.0_wp*ra(l1)*t2a
         vs(i) = vs(i)+t2a*ra(l1)*ra(l1)
      enddo
      ! ---
   enddo

   !      call prmat(6,vs,nat,0,'vs')

end subroutine setvsdq
