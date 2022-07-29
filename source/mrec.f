      subroutine mrec(molcount,xyz,nat,at,molvec)
! molcount: number of total fragments (increased during search)
! xyz: overall Cart. coordinates
! nat: overall number of atoms
! at: atomic number array
! molvec: assignment vector of atom to fragment

      implicit none
      real*8 xyz(3,nat),cn(nat),bond(nat,nat)
      integer nat,molvec(nat),i,molcount,at(nat)
      logical taken(nat)
      molvec=0
      molcount=1
      taken=.false.
      cn=0.0d0
      bond=0.0d0
      call xcoord(nat,at,xyz,cn,bond)
      do i=1,nat
       if(.not.taken(i)) then
         molvec(i)=molcount
         taken(i)=.true.
         call neighbours(i,xyz,at,taken,nat,cn,bond,molvec,molcount)
         molcount=molcount+1
      endif
      enddo
      molcount=molcount-1
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      recursive subroutine neighbours(i,xyz,iat,taken,nat,cn,bond,
     .                                molvec,molcnt)
      implicit none
      real*8 xyz(3,nat),r,tr,xi(3),cntmp,cn(nat),bond(nat,nat)
      integer i,nat, molcnt,molvec(nat),j,iat(nat),icn,k
      logical taken(nat)
      tr=2.0d0
      
      xi(1:3)=xyz(1:3,i) 
      icn=nint(cn(i))
      do k=1,icn
         j=maxloc(bond(:,i),1)
         bond(j,i)=0.0d0
         if (i .eq. j) cycle
         if (.not.taken(j)) then
            molvec(j)=molcnt
            taken(j)=.true.
            call neighbours(j,xyz,iat,taken,nat,cn,bond,molvec,molcnt)
         endif
      enddo
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C compute coordination numbers by adding an inverse damping function
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine xcoord(natoms,iz,xyz,cn,bond)
      implicit none
      integer iz(natoms),natoms,i,j,k1
      real*8 xyz(3,natoms),cn(natoms)
      real*8 cn_thr,bond(natoms,natoms)
      integer iat
      real*8 dx,dy,dz,r,damp,xn,rr,rco,r2,rcovi,rcovj,rcov(94),ff
c D3 values
      data rcov/
     . 0.80628308, 1.15903197, 3.02356173, 2.36845659, 1.94011865,
     . 1.88972601, 1.78894056, 1.58736983, 1.61256616, 1.68815527,
     . 3.52748848, 3.14954334, 2.84718717, 2.62041997, 2.77159820,
     . 2.57002732, 2.49443835, 2.41884923, 4.43455700, 3.88023730,
     . 3.35111422, 3.07395437, 3.04875805, 2.77159820, 2.69600923,
     . 2.62041997, 2.51963467, 2.49443835, 2.54483100, 2.74640188,
     . 2.82199085, 2.74640188, 2.89757982, 2.77159820, 2.87238349,
     . 2.94797246, 4.76210950, 4.20778980, 3.70386304, 3.50229216,
     . 3.32591790, 3.12434702, 2.89757982, 2.84718717, 2.84718717,
     . 2.72120556, 2.89757982, 3.09915070, 3.22513231, 3.17473967,
     . 3.17473967, 3.09915070, 3.32591790, 3.30072128, 5.26603625,
     . 4.43455700, 4.08180818, 3.70386304, 3.98102289, 3.95582657,
     . 3.93062995, 3.90543362, 3.80464833, 3.82984466, 3.80464833,
     . 3.77945201, 3.75425569, 3.75425569, 3.72905937, 3.85504098,
     . 3.67866672, 3.45189952, 3.30072128, 3.09915070, 2.97316878,
     . 2.92277614, 2.79679452, 2.82199085, 2.84718717, 3.32591790,
     . 3.27552496, 3.27552496, 3.42670319, 3.30072128, 3.47709584,
     . 3.57788113, 5.06446567, 4.56053862, 4.20778980, 3.98102289,
     . 3.82984466, 3.85504098, 3.88023730, 3.90543362 /

      cn_thr=400.0d0
      k1=16
      bond=0.0d0
      cn=0.0d0
      do i=1,natoms
        xn=0.0d0
        rcovi=rcov(iz(i))
        do iat=1,natoms
         if(iat.ne.i)then
            dx=xyz(1,iat)-xyz(1,i)
            dy=xyz(2,iat)-xyz(2,i)
            dz=xyz(3,iat)-xyz(3,i)
            r2=dx*dx+dy*dy+dz*dz
            r=sqrt(r2)
            if (r2.gt.cn_thr) cycle
            rcovj=rcov(iz(iat))
c covalent distance in Bohr
            ff=1.00
!           if(iz(i).eq.1.and.iz(iat).eq.1) ff=1.1
            rco=(rcovi+rcovj)*ff   ! this scaling reduces the size of the clusters
            rr=rco/r
c counting function exponential has a better long-range behavior than MHGs inverse damping
            damp=1.d0/(1.d0+exp(-k1*(rr-1.0d0)))
            bond(iat,i)=damp
            xn=xn+damp
         endif
        enddo
        cn(i)=xn
      enddo

      end subroutine xcoord

