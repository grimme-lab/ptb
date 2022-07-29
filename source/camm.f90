
!------------------------------------------------------------------------------------------------------
! distributed atomic multipole moment interactions: all interactions up to r**-3
! energy evaluation
! nat         : # of atoms
! xyz(3,nat)  : cartesian coordinates
! q(nat)      : atomic partial charges
! dipm(3,nat) : cumulative atomic dipole moments (x,y,z)
! qp(6,nat)   : traceless(!) cumulative atomic quadrupole moments (xx,xy,yy,xz,yz,zz)
! gab3,gab5   : damped R**-3 and R**-5 Coulomb laws, dimension: nat*(nat+1)/2
!               multiplication with numerator then leads to R**-2 and R**-3 decay, respectively
! e           : E_AES
!------------------------------------------------------------------------------------------------------

      subroutine aniso_es(nat,at,rab,xyz,q,dipm,qp,e)
      use com
      use parcom
      implicit none          
      integer nat,at(nat)
      real*8 rab(nat*(nat+1)/2)
      real*8 xyz(3,nat),q(nat)
      real*8 dipm(3,nat),qp(6,nat)
      real*8 e           

      real*8 qp1(6),rr(3),dp1(3),rij(3)
      real*8 edd,e01,e02,e11,r0,r5,r3,r2,tt,tt3,q1,qs2,xk,damp
      real*8 ed,eq
      real*8 gab3(nat*(nat+1)/2),gab5(nat*(nat+1)/2)
      integer i,j,k,l,m,ki,kj,kl,lin

      k=0
      do i=1, nat
         do j=1,i-1
            k = k + 1
            r3= rab(k)**3
            r5= r3*rab(k)*rab(k)
            r0= glob_par(19)*(rcov(at(i))+rcov(at(j)))
            damp= 1d0 - 0.5d0*(erf(-1.00*(rab(k)-r0)/r0)+1d0)  
            gab3(k) = damp/r3
            gab5(k) = damp/r5
         enddo
         k=k+1
      enddo

      e=0.0d0
      e01=0.0d0
      e02=0.0d0
      e11=0.0d0
      do i=2,nat 
         q1=q(i)! *glob_par(9)        ! charge scaling
         rr(1:3)=xyz(1:3,i)
         dp1(1:3)=dipm(1:3,i) 
         qp1(1:6)=qp(1:6,i)   
         do j=1,i-1      
            kj=lin(j,i)
            qs2=1d0 !glob_par(9)         ! charge scaling
            rij(1:3)=xyz(1:3,j)-rr(1:3)
            r2=rab(kj)            
            ed=0.0d0
            eq=0.0d0
            edd=0.0d0
!           dipole - charge
            do k=1,3
               ed=ed+q(j)*dp1(k)*rij(k)*qs2
               ed=ed-dipm(k,j)*q1*rij(k)
!              dip-dip & charge-qpole
               do l=1,3
                  kl=lin(l,k)
                  tt=rij(l)*rij(k)
                  tt3=3.0d0*tt
                  eq=eq+q(j)*qp1(kl)*tt*qs2
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
! nat              : # of atoms
! nao              : # of spherical AOs (SAOs)
! s                : overlap matrix
! xyz(3,nat)       : cartesian coordinates
!  pint            : dipole integral matrix, dimension 1-3,nao*(nao+1)/2
!                  : quadrupole integral matrix, dimension 4-9,nao*(nao+1)/2
! s                : exact overlap matrix
! p                : density matrix
! dipm(3,nat)      : cumulative atomic dipole moments (x,y,z)
! qp(6,nat)        : traceless(!) cumulative atomic quadrupole moments (xx,xy,yy,xz,yz,zz)
!------------------------------------------------------------------------------------------------------

subroutine camm(nat,nao,xyz,p,s,pint,dipm,qp)
      use bascom
      implicit none
      integer, intent(in) :: nao,nat
      real*8, intent(in) :: pint(9,nao*(nao+1)/2)
      real*8, intent(in) :: s(nao*(nao+1)/2),p(nao*(nao+1)/2)
      real*8, intent(in) :: xyz(3,nat)
      real*8, intent(out):: dipm(3,nat),qp(6,nat)

      real*8 xk1,xl1,xk2,xl2,pij,tii,tjj
      real*8 pqm,pdmk,pdml,ps,ra(3)

      integer i,j,k,l,ii,jj,ij,kl,kj,lin

      dipm=0.0d0
      qp=0.0d0

      ij=0
      do i=1,nao
         ii=aoat(i)
         ra(1:3)=xyz(1:3,ii)
         do j=1,i-1
            ij=ij+1
            jj=aoat(j)
            pij=p(ij)
            ps=pij*s(ij)
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
         ij=ij+1
         pij=p(ij)
         ps=pij*s(ij)
         kl=0
         !  the qpint is stored as xx,yy,zz,xy,xz,yz (from integral routine)
         !  when doing the Mulliken population, we switch to lin-compatible sorting
         !  i,e. xx,xy,yy,xz,yz,zz 
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
      enddo 

end subroutine camm   

