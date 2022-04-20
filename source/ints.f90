
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sint(nat,nao,at,xyz,rab,s,norm)
      use bascom
      implicit none          
      integer, intent(in)  :: nao,nat,at(nat)
      real*8, intent(in)   :: xyz(3,nat)
      real*8, intent(in)   :: rab(nat*(nat+1)/2)
      real*8, intent(out)  :: s(nao*(nao+1)/2)
      real*8, intent(out)  :: norm(nao)

      real*8 ss(6,6)
      real*8 tmp1,tmp2,tmp3,tmp4,intcut,s00,apb
      real*8 r1,r2,t1,t2,t3,t4,f,ci,cc,cj,alpi,rab2,ab,est
      real*8 rc(3),ra(3),rb(3),p(3),g(3),alpj,va(1),alpc,ccc
      integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kl,km,mi,mj,ij,lin,jshmax,ksh
      integer ip,jp,iat,jat,iatyp,jatyp,ish,jsh,icao,jcao,iao,jao,isao
      integer ishtyp,jshtyp,naoi,naoj,li,lj,iprim,jprim,iijj,iptyp,jptyp
      integer llao(0:3),llao2(0:3),itt(0:3),mli,mlj,nuc,abcnt
      data llao /1,3,6,10/
      data llao2/1,3,5,7 /
      data itt  /0,1,4,10/
      real*8, parameter :: PI = 2.0d0/3.1415926535897932384626433832795029d0
      real*8, parameter :: PI2= 2.0d0/PI

      intcut=30. ! cut-off
      
      s =0.0d0
      do iat=1,nat
         ra(1:3)=xyz(1:3,iat)
         iatyp=at(iat)
!        write (*,*), 'i',iat, iatyp, bas_nsh(iatyp)
         do jat=1,iat
            abcnt=lin(iat,jat)
            rb(1:3)=xyz(1:3,jat)
            jatyp=at(jat)
            rab2=rab(abcnt)**2               
!c          ints < 1.d-9 for RAB > 40 Bohr            
            if(rab2.gt.1600) cycle
            do ish=1,bas_nsh(iatyp)
               ishtyp=bas_lsh(ish,iatyp)
               icao=caoshell(ish,iat)
               naoi=llao(ishtyp)
               iptyp=itt(ishtyp)
               jshmax=bas_nsh(jatyp)
!              jshells
               do jsh=1,jshmax
                  jshtyp=bas_lsh(jsh,jatyp)             
                  jcao=caoshell(jsh,jat)
                  naoj=llao(jshtyp)
                  jptyp=itt(jshtyp)
                  ! we go through the primitives (beacause the screening is the same for all of them) 
                  ss=0.0d0
                  do ip=1,prim_npr(icao+1)
                     iprim=ip+prim_count(icao+1)
                     alpi=prim_exp(iprim) ! exponent the same for each l component
                     do jp=1,prim_npr(jcao+1)
                        jprim=jp+prim_count(jcao+1)
                        alpj=prim_exp(jprim) ! exponent the same for each l component
                        apb = alpi + alpj 
                        ab  =1.0d0/apb        
                        est =rab2*alpi*alpj*ab
                        if(est.gt.intcut) cycle

!                       now compute integrals for different components of i(e.g., px,py,pz)
! S and T
                        do mli=1,naoi
                           iprim=ip+prim_count(icao+mli)
                           ci=prim_cnt(iprim) ! coefficients NOT the same (contain CAO2SAO lin. comb. coefficients)
                           do mlj=1,naoj
                             jprim=jp+prim_count(jcao+mlj)
                             cc=prim_cnt(jprim)*ci
                             call propa_s  (nat,xyz,ra,rb,alpi,alpj,iptyp+mli,jptyp+mlj,va)
                             ss(mlj,mli)=ss(mlj,mli)+va(1)*cc
                           enddo ! mlj
                        enddo ! mli
                     enddo ! jp
                  enddo ! ip 

                  !transform from CAO to SAO
                  call dtrf2(ss,ishtyp,jshtyp)
                  do ii=1,llao2(ishtyp)
                     iao=ii+aoshell(ish,iat)
                     do jj=1,llao2(jshtyp)
                        jao=jj+aoshell(jsh,jat)
                        ij = lin(jao,iao)
                        s(ij)=ss(jj,ii)
                     enddo
                  enddo
               enddo!jsh
            enddo!ish
         enddo!jat
      enddo!iat

      ij = 0
      do i=1,nao
         ij = ij + i
         norm(i)=1./sqrt(s(ij))
      enddo

      ij = 0
      do i=1,nao
         do j=1,i   
            ij = ij + 1
            tmp1=norm(i)*norm(j)
            s(ij)=s(ij)*tmp1
         enddo
      enddo

end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine stint(nat,nao,at,xyz,rab,s,t,norm)
      use bascom
      implicit none          
      integer, intent(in)  :: nao,nat,at(nat)
      real*8, intent(in)   :: xyz(3,nat)
      real*8, intent(in)   :: rab(nat*(nat+1)/2)
      real*8, intent(out)  :: s(nao*(nao+1)/2)
      real*8, intent(out)  :: t(nao*(nao+1)/2)
      real*8, intent(out)  :: norm(nao)

      real*8 ss(6,6), tt(6,6), va2(6,6)
      real*8 tmp1,tmp2,tmp3,tmp4,intcut,s00,apb
      real*8 r1,r2,t1,t2,t3,t4,f,ci,cc,cj,alpi,rab2,ab,est
      real*8 rc(3),ra(3),rb(3),p(3),g(3),alpj,va(2+nat),alpc,ccc
      integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kl,km,mi,mj,ij,lin,jshmax,ksh
      integer ip,jp,iat,jat,iatyp,jatyp,ish,jsh,icao,jcao,iao,jao,isao
      integer ishtyp,jshtyp,naoi,naoj,li,lj,iprim,jprim,iijj,iptyp,jptyp
      integer llao(0:3),llao2(0:3),itt(0:3),mli,mlj,nuc,abcnt
      data llao /1,3,6,10/
      data llao2/1,3,5,7 /
      data itt  /0,1,4,10/
      real*8, parameter :: PI = 2.0d0/3.1415926535897932384626433832795029d0
      real*8, parameter :: PI2= 2.0d0/PI

      intcut=30. ! agressive cut-off
      
      s =0.0d0
      t =0.0d0
      do iat=1,nat
         ra(1:3)=xyz(1:3,iat)
         iatyp=at(iat)
         do jat=1,iat
            abcnt=lin(iat,jat)
            rb(1:3)=xyz(1:3,jat)
            jatyp=at(jat)
            rab2=rab(abcnt)**2               
!c          ints < 1.d-9 for RAB > 40 Bohr            
            if(rab2.gt.900) cycle
            do ish=1,bas_nsh(iatyp)
               ishtyp=bas_lsh(ish,iatyp)
               icao=caoshell(ish,iat)
               naoi=llao(ishtyp)
               iptyp=itt(ishtyp)
               jshmax=bas_nsh(jatyp)
!              jshells
               do jsh=1,jshmax
                  jshtyp=bas_lsh(jsh,jatyp)             
                  jcao=caoshell(jsh,jat)
                  naoj=llao(jshtyp)
                  jptyp=itt(jshtyp)
                  ! we go through the primitives (beacause the screening is the same for all of them) 
                  ss=0.0d0
                  tt=0.0d0
                  do ip=1,prim_npr(icao+1)
                     iprim=ip+prim_count(icao+1)
                     alpi=prim_exp(iprim) ! exponent the same for each l component
                     do jp=1,prim_npr(jcao+1)
                        jprim=jp+prim_count(jcao+1)
                        alpj=prim_exp(jprim) ! exponent the same for each l component
                        apb = alpi + alpj 
                        ab  =1.0d0/apb        
                        est =rab2*alpi*alpj*ab
                        if(est.gt.intcut) cycle

!                       now compute integrals for different components of i(e.g., px,py,pz)
! S and T
                        do mli=1,naoi
                           iprim=ip+prim_count(icao+mli)
                           ci=prim_cnt(iprim) ! coefficients NOT the same (contain CAO2SAO lin. comb. coefficients)
                           do mlj=1,naoj
                             jprim=jp+prim_count(jcao+mlj)
                             cc=prim_cnt(jprim)*ci
                             call propa_st (nat,xyz,ra,rb,alpi,alpj,iptyp+mli,jptyp+mlj,va)
                             ss(mlj,mli)=ss(mlj,mli)+va(1)*cc
                             tt(mlj,mli)=tt(mlj,mli)+va(2)*cc
                           enddo ! mlj
                        enddo ! mli
                     enddo ! jp
                  enddo ! ip 

                  !transform from CAO to SAO
                  call dtrf2(ss,ishtyp,jshtyp)
                  call dtrf2(tt,ishtyp,jshtyp)
                  do ii=1,llao2(ishtyp)
                     iao=ii+aoshell(ish,iat)
                     do jj=1,llao2(jshtyp)
                        jao=jj+aoshell(jsh,jat)
                        ij = lin(jao,iao)
                        s(ij)=ss(jj,ii)
                        t(ij)=tt(jj,ii)
                     enddo
                  enddo
               enddo!jsh
            enddo!ish
         enddo!jat
      enddo!iat

      ij = 0
      do i=1,nao
         ij = ij + i
         norm(i)=1./sqrt(s(ij))
      enddo

      ij = 0
      do i=1,nao
         do j=1,i   
            ij = ij + 1
            tmp1=norm(i)*norm(j)
            s(ij)=s(ij)*tmp1
            t(ij)=t(ij)*tmp1
         enddo
      enddo

end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dipmom(nat,nao,at,xyz,z,rab,norm,p,pnt,dip)
      use bascom
      implicit none          
      integer, intent(in)  :: nao,nat,at(nat)
      real*8, intent(in)   :: xyz(3,nat)
      real*8, intent(in)   :: z(nat)
      real*8, intent(in)   :: rab(nat*(nat+1)/2)
      real*8, intent(in)   :: norm(nao)
      real*8, intent(in)   :: p(nao*(nao+1)/2)
      real*8, intent(in)   :: pnt(3) ! reference point
      real*8, intent(out)  :: dip(3)

      real*8 ss3(6,6,3)
      real*8 tmp1,tmp2,tmp3,tmp4,intcut,s00,apb
      real*8 r1,r2,t1,t2,t3,t4,f,ci,cc,cj,alpi,rab2,ab,est
      real*8 rc(3),ra(3),rb(3),alpj,va(3),alpc,ccc
      integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kl,km,mi,mj,ij,lin,jshmax,ksh
      integer ip,jp,iat,jat,iatyp,jatyp,ish,jsh,icao,jcao,iao,jao,isao
      integer ishtyp,jshtyp,iptyp,jptyp,naoi,naoj,li,lj,iprim,jprim,iijj
      integer llao(0:3),llao2(0:3),itt(0:3),mli,mlj,nuc,abcnt
      data llao /1,3,6,10/
      data llao2/1,3,5,7 /
      data itt  /0,1,4,10/
      real*8,allocatable :: d(:,:) 

      allocate(d(3,nao*(nao+1)/2))
      intcut=25. 
      
      d =0.0d0
      do iat=1,nat
         ra(1:3)=xyz(1:3,iat)
         iatyp=at(iat)
!        write (*,*), 'i',iat, iatyp, bas_nsh(iatyp)
         do jat=1,iat
            abcnt=lin(iat,jat)
            rb(1:3)=xyz(1:3,jat)
            jatyp=at(jat)
            rab2=rab(abcnt)**2               
!c          ints < 1.d-9 for RAB > 40 Bohr            
            if(rab2.gt.900) cycle
            do ish=1,bas_nsh(iatyp)
               ishtyp=bas_lsh(ish,iatyp)
               icao=caoshell(ish,iat)
               naoi=llao(ishtyp)
               iptyp=itt(ishtyp)
               jshmax=bas_nsh(jatyp)
!              jshells
               do jsh=1,jshmax
                  jshtyp=bas_lsh(jsh,jatyp)             
                  jcao=caoshell(jsh,jat)
                  naoj=llao(jshtyp)
                  jptyp=itt(jshtyp)
                  ! we go through the primitives (beacause the screening is the same for all of them) 
                  ss3=0.0d0
                  do ip=1,prim_npr(icao+1)
                     iprim=ip+prim_count(icao+1)
                     alpi=prim_exp(iprim) ! exponent the same for each l component
                     do jp=1,prim_npr(jcao+1)
                        jprim=jp+prim_count(jcao+1)
                        alpj=prim_exp(jprim) ! exponent the same for each l component
                        apb = alpi + alpj 
                        ab  =1.0d0/apb        
                        est =rab2*alpi*alpj*ab
                        if(est.gt.intcut) cycle
!                       now compute integrals for different components of i(e.g., px,py,pz)
                        do mli=1,naoi
                           iprim=ip+prim_count(icao+mli)
                           ci=prim_cnt(iprim) ! coefficients NOT the same (contain CAO2SAO lin. comb. coefficients)
                           do mlj=1,naoj
                             jprim=jp+prim_count(jcao+mlj)
                             cc=prim_cnt(jprim)*ci
                             call propa_dip(pnt,ra,rb,alpi,alpj,iptyp+mli,jptyp+mlj,va)
                             do ksh=1,3
                                ss3(mlj,mli,ksh)=ss3(mlj,mli,ksh)+va(ksh)*cc
                             enddo
                           enddo ! mlj
                        enddo ! mli
                     enddo ! jp
                  enddo ! ip 
                  !transform from CAO to SAO
                  do ksh=1,3     
                     call dtrf2(ss3(1,1,ksh),ishtyp,jshtyp)
                  enddo
                  do ii=1,llao2(ishtyp)
                     iao=ii+aoshell(ish,iat)
                     do jj=1,llao2(jshtyp)
                        jao=jj+aoshell(jsh,jat)
                        ij = lin(jao,iao)
                        d(1:3,ij)=ss3(jj,ii,1:3)
                     enddo
                  enddo
               enddo!jsh
            enddo!ish
         enddo!jat
      enddo!iat

      va = 0

      ij = 0
      do i=1,nao
         do j=1,i-1
            ij = ij + 1
            tmp1=norm(i)*norm(j)
            va(1:3)=va(1:3)+d(1:3,ij)*tmp1*P(ij)
         enddo
         ij = ij + 1
         va(1:3)=va(1:3)+d(1:3,ij)*norm(i)*norm(i)*P(ij)*0.50d0
      enddo
      va = va * 2.d0

      do i=1,nat
         va(1:3)=va(1:3)+(xyz(1:3,i)-pnt(1:3))*z(i)
      enddo

      dip = va

end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine stvint(nat,nao,at,xyz,rab,z,s,t,v,norm)
      use bascom
      implicit none          
      integer, intent(in)  :: nao,nat,at(nat)
      real*8, intent(in)   :: xyz(3,nat)
      real*8, intent(in)   :: rab(nat*(nat+1)/2)
      real*8, intent(in)   :: z(nat)
      real*8, intent(out)  :: s(nao*(nao+1)/2)
      real*8, intent(out)  :: t(nao*(nao+1)/2)
      real*8, intent(out)  :: v(nao*(nao+1)/2)
      real*8, intent(out)  :: norm(nao)

      real*8 ss(6,6), tt(6,6), vv(6,6), va2(6,6)
      real*8 tmp1,tmp2,tmp3,tmp4,intcut,s00,apb
      real*8 r1,r2,t1,t2,t3,t4,f,ci,cc,cj,alpi,rab2,ab,est
      real*8 rc(3),ra(3),rb(3),p(3),g(3),alpj,va(2+nat),alpc,ccc
      real*8 abc,eabc,rpq2,s000
      integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kl,km,mi,mj,ij,lin,jshmax,ksh
      integer ip,jp,iat,jat,iatyp,jatyp,ish,jsh,icao,jcao,iao,jao,isao
      integer ishtyp,jshtyp,iptyp,jptyp,naoi,naoj,li,lj,iprim,jprim,iijj
      integer llao(0:3),llao2(0:3),itt(0:3),mli,mlj,nuc,abcnt
      data llao /1,3,6,10/
      data llao2/1,3,5,7 /
      data itt  /0,1,4,10/
      real*8, parameter :: PI = 2.0d0/3.1415926535897932384626433832795029d0
      real*8, parameter :: PI2= 2.0d0/PI

      intcut=25. ! agressive cut-off
      
      s =0.0d0
      t =0.0d0
      v =0.0d0
      do iat=1,nat
         ra(1:3)=xyz(1:3,iat)
         iatyp=at(iat)
!        write (*,*), 'i',iat, iatyp, bas_nsh(iatyp)
         do jat=1,iat
            abcnt=lin(iat,jat)
            rb(1:3)=xyz(1:3,jat)
            jatyp=at(jat)
            rab2=rab(abcnt)**2               
!c          ints < 1.d-9 for RAB > 40 Bohr            
            if(rab2.gt.2500) cycle
!           write(*,*) 'j',jat,jatyp,bas_nsh(jatyp)
            do ish=1,bas_nsh(iatyp)
               ishtyp=bas_lsh(ish,iatyp)
               icao=caoshell(ish,iat)
               naoi=llao(ishtyp)
               iptyp=itt(ishtyp)
!              write(*,*) 'i',iat,iatyp,ish,ishtyp,icao,naoi
               iptyp=itt(ishtyp)
               jshmax=bas_nsh(jatyp)
!              jshells
               do jsh=1,jshmax
                  jshtyp=bas_lsh(jsh,jatyp)             
                  jcao=caoshell(jsh,jat)
                  naoj=llao(jshtyp)
                  jptyp=itt(jshtyp)
!                 write(*,*) 'j',jat,jatyp,jsh,jshtyp,jcao,naoj
                  ! we go through the primitives (beacause the screening is the same for all of them) 
                  ss=0.0d0
                  tt=0.0d0
                  vv=0.0d0
                  do ip=1,prim_npr(icao+1)
                     iprim=ip+prim_count(icao+1)
                     alpi=prim_exp(iprim) ! exponent the same for each l component
                     do jp=1,prim_npr(jcao+1)
                        jprim=jp+prim_count(jcao+1)
                        alpj=prim_exp(jprim) ! exponent the same for each l component
                        apb = alpi + alpj 
                        ab  =1.0d0/apb        
                        est =rab2*alpi*alpj*ab
                        if(est.gt.intcut) cycle

!                       now compute integrals for different components of i(e.g., px,py,pz)
! S and T
                        do mli=1,naoi
                           iprim=ip+prim_count(icao+mli)
                           ci=prim_cnt(iprim) ! coefficients NOT the same (contain CAO2SAO lin. comb. coefficients)
                           do mlj=1,naoj
                             jprim=jp+prim_count(jcao+mlj)
                             cc=prim_cnt(jprim)*ci
!                            call propa_st (nat,xyz,ra,rb,alpi,alpj,iptyp+mli,jptyp+mlj,va)
                             call propa_stv(nat,xyz,ra,rb,alpi,alpj,iptyp+mli,jptyp+mlj,va)
                             ss(mlj,mli)=ss(mlj,mli)+va(1)*cc
                             tt(mlj,mli)=tt(mlj,mli)+va(2)*cc
                             do nuc=1,nat
                             vv(mlj,mli)=vv(mlj,mli)+va(nuc+2)*cc*z(nuc)
                             enddo
                           enddo ! mlj
                        enddo ! mli
                     enddo ! jp
                  enddo ! ip 

                  !transform from CAO to SAO
                  call dtrf2(ss,ishtyp,jshtyp)
                  call dtrf2(tt,ishtyp,jshtyp)
                  call dtrf2(vv,ishtyp,jshtyp)
                  do ii=1,llao2(ishtyp)
                     iao=ii+aoshell(ish,iat)
                     do jj=1,llao2(jshtyp)
                        jao=jj+aoshell(jsh,jat)
                        ij = lin(jao,iao)
                        s(ij)=ss(jj,ii)
                        t(ij)=tt(jj,ii)
                        v(ij)=vv(jj,ii)
                     enddo
                  enddo
               enddo!jsh
            enddo!ish
         enddo!jat
      enddo!iat

      ij = 0
      do i=1,nao
         ij = ij + i
         norm(i)=1./sqrt(s(ij))
      enddo

      ij = 0
      do i=1,nao
         do j=1,i   
            ij = ij + 1
            tmp1=norm(i)*norm(j)
            s(ij)=s(ij)*tmp1
            t(ij)=t(ij)*tmp1
            v(ij)=v(ij)*tmp1
         enddo
      enddo

end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dipint(nat,nao,at,xyz,rab,norm,pnt,d)
      use bascom
      implicit none          
      integer, intent(in)  :: nao,nat,at(nat)
      real*8, intent(in)   :: xyz(3,nat)
      real*8, intent(in)   :: rab(nat*(nat+1)/2)
      real*8, intent(in)   :: norm(nao)
      real*8, intent(in)   :: pnt(3) ! reference point
      real*8, intent(out)  :: d(nao*(nao+1)/2,3)

      real*8 ss3(6,6,3)
      real*8 tmp1,tmp2,tmp3,tmp4,intcut,s00,apb
      real*8 r1,r2,t1,t2,t3,t4,f,ci,cc,cj,alpi,rab2,ab,est
      real*8 rc(3),ra(3),rb(3),alpj,va(3),alpc,ccc
      integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kl,km,mi,mj,ij,lin,jshmax,ksh
      integer ip,jp,iat,jat,iatyp,jatyp,ish,jsh,icao,jcao,iao,jao,isao
      integer ishtyp,jshtyp,iptyp,jptyp,naoi,naoj,li,lj,iprim,jprim,iijj
      integer llao(0:3),llao2(0:3),itt(0:3),mli,mlj,nuc,abcnt
      data llao /1,3,6,10/
      data llao2/1,3,5,7 /
      data itt  /0,1,4,10/

      intcut=20. 
      
      d =0.0d0
      do iat=1,nat
         ra(1:3)=xyz(1:3,iat)
         iatyp=at(iat)
!        write (*,*), 'i',iat, iatyp, bas_nsh(iatyp)
         do jat=1,iat
            abcnt=lin(iat,jat)
            rb(1:3)=xyz(1:3,jat)
            jatyp=at(jat)
            rab2=rab(abcnt)**2               
!c          ints < 1.d-9 for RAB > 40 Bohr            
            if(rab2.gt.400) cycle
            do ish=1,bas_nsh(iatyp)
               ishtyp=bas_lsh(ish,iatyp)
               icao=caoshell(ish,iat)
               naoi=llao(ishtyp)
               iptyp=itt(ishtyp)
               jshmax=bas_nsh(jatyp)
!              jshells
               do jsh=1,jshmax
                  jshtyp=bas_lsh(jsh,jatyp)             
                  jcao=caoshell(jsh,jat)
                  naoj=llao(jshtyp)
                  jptyp=itt(jshtyp)
                  ! we go through the primitives (beacause the screening is the same for all of them) 
                  ss3=0.0d0
                  do ip=1,prim_npr(icao+1)
                     iprim=ip+prim_count(icao+1)
                     alpi=prim_exp(iprim) ! exponent the same for each l component
                     do jp=1,prim_npr(jcao+1)
                        jprim=jp+prim_count(jcao+1)
                        alpj=prim_exp(jprim) ! exponent the same for each l component
                        apb = alpi + alpj 
                        ab  =1.0d0/apb        
                        est =rab2*alpi*alpj*ab
                        if(est.gt.intcut) cycle
!                       now compute integrals for different components of i(e.g., px,py,pz)
                        do mli=1,naoi
                           iprim=ip+prim_count(icao+mli)
                           ci=prim_cnt(iprim) ! coefficients NOT the same (contain CAO2SAO lin. comb. coefficients)
                           do mlj=1,naoj
                             jprim=jp+prim_count(jcao+mlj)
                             cc=prim_cnt(jprim)*ci
                             call propa_dip(pnt,ra,rb,alpi,alpj,iptyp+mli,jptyp+mlj,va)
                             do ksh=1,3
                                ss3(mlj,mli,ksh)=ss3(mlj,mli,ksh)+va(ksh)*cc
                             enddo
                           enddo ! mlj
                        enddo ! mli
                     enddo ! jp
                  enddo ! ip 
                  !transform from CAO to SAO
                  do ksh=1,3     
                     call dtrf2(ss3(1,1,ksh),ishtyp,jshtyp)
                  enddo
                  do ii=1,llao2(ishtyp)
                     iao=ii+aoshell(ish,iat)
                     do jj=1,llao2(jshtyp)
                        jao=jj+aoshell(jsh,jat)
                        ij = lin(jao,iao)
                        d(ij,1:3)=ss3(jj,ii,1:3)
                     enddo
                  enddo
               enddo!jsh
            enddo!ish
         enddo!jat
      enddo!iat

      do k=1,3
      ij = 0
      do i=1,nao
         do j=1,i
            ij = ij + 1
            tmp1=norm(i)*norm(j)
            d(ij,k)=d(ij,k)*tmp1
         enddo
      enddo
      enddo

end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! second moment, x^2,y^2,z^2 components only
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine secint(nat,nao,at,xyz,rab,norm,pnt,d)
      use bascom
      implicit none          
      integer, intent(in)  :: nao,nat,at(nat)
      real*8, intent(in)   :: xyz(3,nat)
      real*8, intent(in)   :: rab(nat*(nat+1)/2)
      real*8, intent(in)   :: norm(nao)
      real*8, intent(in)   :: pnt(3) ! reference point
      real*8, intent(out)  :: d(nao*(nao+1)/2,3)

      real*8 ss3(6,6,3)
      real*8 tmp1,tmp2,tmp3,tmp4,intcut,s00,apb
      real*8 r1,r2,t1,t2,t3,t4,f,ci,cc,cj,alpi,rab2,ab,est
      real*8 rc(3),ra(3),rb(3),alpj,va(6),alpc,ccc
      integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kl,km,mi,mj,ij,lin,jshmax,ksh
      integer ip,jp,iat,jat,iatyp,jatyp,ish,jsh,icao,jcao,iao,jao,isao
      integer ishtyp,jshtyp,iptyp,jptyp,naoi,naoj,li,lj,iprim,jprim,iijj
      integer llao(0:3),llao2(0:3),itt(0:3),mli,mlj,nuc,abcnt
      data llao /1,3,6,10/
      data llao2/1,3,5,7 /
      data itt  /0,1,4,10/

      intcut=20. 
      
      d =0.0d0
      do iat=1,nat
         ra(1:3)=xyz(1:3,iat)
         iatyp=at(iat)
!        write (*,*), 'i',iat, iatyp, bas_nsh(iatyp)
         do jat=1,iat
            abcnt=lin(iat,jat)
            rb(1:3)=xyz(1:3,jat)
            jatyp=at(jat)
            rab2=rab(abcnt)**2               
!c          ints < 1.d-9 for RAB > 40 Bohr            
            if(rab2.gt.400) cycle
            do ish=1,bas_nsh(iatyp)
               ishtyp=bas_lsh(ish,iatyp)
               icao=caoshell(ish,iat)
               naoi=llao(ishtyp)
               iptyp=itt(ishtyp)
               jshmax=bas_nsh(jatyp)
!              jshells
               do jsh=1,jshmax
                  jshtyp=bas_lsh(jsh,jatyp)             
                  jcao=caoshell(jsh,jat)
                  naoj=llao(jshtyp)
                  jptyp=itt(jshtyp)
                  ! we go through the primitives (beacause the screening is the same for all of them) 
                  ss3=0.0d0
                  do ip=1,prim_npr(icao+1)
                     iprim=ip+prim_count(icao+1)
                     alpi=prim_exp(iprim) ! exponent the same for each l component
                     do jp=1,prim_npr(jcao+1)
                        jprim=jp+prim_count(jcao+1)
                        alpj=prim_exp(jprim) ! exponent the same for each l component
                        apb = alpi + alpj 
                        ab  =1.0d0/apb        
                        est =rab2*alpi*alpj*ab
                        if(est.gt.intcut) cycle
!                       now compute integrals for different components of i(e.g., px,py,pz)
                        do mli=1,naoi
                           iprim=ip+prim_count(icao+mli)
                           ci=prim_cnt(iprim) ! coefficients NOT the same (contain CAO2SAO lin. comb. coefficients)
                           do mlj=1,naoj
                             jprim=jp+prim_count(jcao+mlj)
                             cc=prim_cnt(jprim)*ci
                             call propa_sec(pnt,ra,rb,alpi,alpj,iptyp+mli,jptyp+mlj,va)
                             do ksh=1,3
                                ss3(mlj,mli,ksh)=ss3(mlj,mli,ksh)+va(ksh)*cc
                             enddo
                           enddo ! mlj
                        enddo ! mli
                     enddo ! jp
                  enddo ! ip 
                  !transform from CAO to SAO
                  do ksh=1,3     
                     call dtrf2(ss3(1,1,ksh),ishtyp,jshtyp)
                  enddo
                  do ii=1,llao2(ishtyp)
                     iao=ii+aoshell(ish,iat)
                     do jj=1,llao2(jshtyp)
                        jao=jj+aoshell(jsh,jat)
                        ij = lin(jao,iao)
                        d(ij,1:3)=ss3(jj,ii,1:3)
                     enddo
                  enddo
               enddo!jsh
            enddo!ish
         enddo!jat
      enddo!iat

      do k=1,3
      ij = 0
      do i=1,nao
         do j=1,i
            ij = ij + 1
            tmp1=norm(i)*norm(j)
            d(ij,k)=d(ij,k)*tmp1
         enddo
      enddo
      enddo

end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! true dipole moment
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dipmom2(nat,nao,xyz,z,norm,p,d,pnt,dip)
      use bascom
      implicit none          
      integer, intent(in)  :: nao,nat
      real*8, intent(in)   :: xyz(3,nat)
      real*8, intent(in)   :: z(nat)
      real*8, intent(in)   :: norm(nao)
      real*8, intent(in)   :: p(nao*(nao+1)/2)
      real*8, intent(in)   :: d(nao*(nao+1)/2,3)
      real*8, intent(in)   :: pnt(3) ! reference point
      real*8, intent(out)  :: dip(3)

      integer i,j,k,ij
      real*8 va, tmp1

      do k=1,3
      va = 0
      ij = 0
      do i=1,nao
         do j=1,i-1
            ij = ij + 1
            va=va+d(ij,k)*P(ij)
         enddo
         ij = ij + 1
         va=va+d(ij,k)*P(ij)*0.50d0
      enddo
      va = va * 2.d0
      do i=1,nat
         va=va+(xyz(k,i)-pnt(k))*z(i)
      enddo
      dip(k) = va
      enddo

end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! just electronic dipole moment
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dipmom3(nat,nao,xyz,z,norm,p,d,dip)
      use bascom
      implicit none          
      integer, intent(in)  :: nao,nat
      real*8, intent(in)   :: xyz(3,nat)
      real*8, intent(in)   :: z(nat)
      real*8, intent(in)   :: norm(nao)
      real*8, intent(in)   :: p(nao*(nao+1)/2)
      real*8, intent(in)   :: d(nao*(nao+1)/2,3)
      real*8, intent(out)  :: dip(3)

      integer i,j,k,ij
      real*8 va, tmp1

      do k=1,3
      va = 0
      ij = 0
      do i=1,nao
         do j=1,i-1
            ij = ij + 1
            va=va+d(ij,k)*P(ij)
         enddo
         ij = ij + 1
         va=va+d(ij,k)*P(ij)*0.50d0
      enddo
      va = va * 2.d0
      dip(k) = -va
      enddo

end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! second moment X^2, Y^2, Z^2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine secmom(nat,nao,xyz,z,norm,p,d,pnt,sec)
      use bascom
      implicit none          
      integer, intent(in)  :: nao,nat
      real*8, intent(in)   :: xyz(3,nat)
      real*8, intent(in)   :: z(nat)
      real*8, intent(in)   :: norm(nao)
      real*8, intent(in)   :: p(nao*(nao+1)/2)
      real*8, intent(in)   :: d(nao*(nao+1)/2,3)
      real*8, intent(in)   :: pnt(3) ! reference point
      real*8, intent(out)  :: sec(3)

      integer i,j,k,ij
      real*8 va, tmp1

      do k=1,3
      va = 0
      ij = 0
      do i=1,nao
         do j=1,i-1
            ij = ij + 1
            va=va+d(ij,k)*P(ij)
         enddo
         ij = ij + 1
         va=va+d(ij,k)*P(ij)*0.50d0
      enddo
      va = va * 2.d0
      do i=1,nat
         va=va-z(i)*(xyz(k,i)-pnt(k))**2
      enddo
      sec(k) = va
      enddo

end
