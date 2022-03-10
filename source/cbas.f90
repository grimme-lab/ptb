
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! simple ECP 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calcvecp(n,nao,at,xyz,rab,norm,v)
      use cbascom
      use  bascom
      use  parcom
      use gtb_la, only : la_gemm
      implicit none          
      integer, intent(in)   :: nao,n,at(n)
      real*8,  intent(in)   :: xyz(3,n)
      real*8,  intent(in)   :: rab(n*(n+1)/2) 
      real*8,  intent(in)   :: norm(nao)    
      real*8,  intent(out)  :: v(nao*(nao+1)/2)    

      integer i,j,k,l,m,nl,nn,atn,jsh,llao2(0:3),ia,ib
      data llao2/1,3,5,7 /
      real*8 vecp, ddot
      real*8,allocatable :: Scv(:,:), stmp(:,:), xtmp(:,:)

      v = 0 

      if(cnsao.eq.0) return

      allocate(Scv(cnsao,nao),stmp(cnsao,nao),xtmp(nao,nao))

      call csint(n,nao,at,xyz,rab,norm,Scv) ! core val overlap ints

!     N^2 step
      do i=1,nao                          
         m=0 
         do nl=1,ncorelist                     ! all atoms with core
            nn=corelist(nl)
            atn=at(nn)
            do jsh=1,cbas_nsh(atn)             ! core shells of atom nn
               do l=1,llao2(cbas_lsh(jsh,atn)) ! AOs of core shell jsh
                  m = m + 1
                  stmp(m,i)=-clev(jsh,atn) * Scv(m,i) * shell_cnf1(10,atn)
               enddo
            enddo
         enddo
      enddo

!     N^3 step
      call la_gemm('T','N',nao,nao,cnsao,1.0d0,Scv,cnsao,stmp,cnsao,0.0d0,xtmp,nao)
      k = 0 
      do i=1, nao 
         do j=1, i
            k = k + 1 
            v(k) = xtmp(j,i) 
         enddo
      enddo

      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! core - valence overlap matrix in SAO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine csint(nat,nao,at,xyz,rab,norm,s)
      use bascom
      use cbascom
      implicit none          
      integer, intent(in)  :: nao,nat,at(nat)
      real*8, intent(in)   :: xyz(3,nat)
      real*8, intent(in)   :: rab(nat*(nat+1)/2)
      real*8, intent(in)   :: norm(nao)        
      real*8, intent(out)  :: s(cnsao,nao)        

      real*8 ss(6,6)
      real*8 tmp1,tmp2,tmp3,tmp4,intcut,s00,apb
      real*8 r1,r2,t1,t2,t3,t4,f,ci,cc,cj,alpi,rab2,ab,est
      real*8 rc(3),ra(3),rb(3),p(3),g(3),alpj,va(1),alpc,ccc
      integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kl,km,mi,mj,ij,lin,jshmax,ksh
      integer ip,jp,iat,jat,iatyp,jatyp,ish,jsh,icao,jcao,iao,jao,isao
      integer ishtyp,jshtyp,naoi,naoj,li,lj,iprim,jprim,iijj,iptyp,jptyp
      integer llao(0:3),llao2(0:3),itt(0:3),mli,mlj,abcnt,nl
      data llao /1,3,6,10/
      data llao2/1,3,5,7 /
      data itt  /0,1,4,10/

      intcut=20. ! cut-off
      
      s =0.0d0
      do iat=1,nat
         ra(1:3)=xyz(1:3,iat)
         iatyp=at(iat)
         do nl=1,ncorelist
            jat=corelist(nl)
            abcnt=lin(iat,jat)
            rb(1:3)=xyz(1:3,jat)
            jatyp=at(jat)
            rab2=rab(abcnt)**2               
!c          ints < 1.d-9 for RAB > 40 Bohr            
            if(rab2.gt.1600) cycle
            do ish=1,bas_nsh(iatyp)              !val
               ishtyp=bas_lsh(ish,iatyp)
               icao=caoshell(ish,iat)
               naoi=llao(ishtyp)
               iptyp=itt(ishtyp)
               do jsh=1,cbas_nsh(jatyp)          !core
                  jshtyp=cbas_lsh(jsh,jatyp)             
                  jcao=caoshellc(jsh,jat)
                  naoj=llao(jshtyp)
                  jptyp=itt(jshtyp)
                  ss=0.0d0 
                  do ip=1,prim_npr(icao+1)       !val
                     iprim=ip+prim_count(icao+1)
                     alpi=prim_exp(iprim) 
                     do jp=1,cprim_npr(jcao+1)   !core
                        jprim=jp+cprim_count(jcao+1)
                        alpj=cprim_exp(jprim) 
                        apb = alpi + alpj 
                        ab  =1.0d0/apb        
                        est =rab2*alpi*alpj*ab
                        if(est.gt.intcut) cycle
                        do mli=1,naoi            !val
                           iprim=ip+prim_count(icao+mli)
                           ci=prim_cnt(iprim)    
                           do mlj=1,naoj         !core
                             jprim=jp+cprim_count(jcao+mlj)
                             cc=cprim_cnt(jprim)*ci
                             call propa_s(nat,xyz,ra,rb,alpi,alpj,iptyp+mli,jptyp+mlj,va)
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
                        jao=jj+aoshellc(jsh,jat)
                        s(jao,iao)=ss(jj,ii) * norm(iao) ! core STOs are properly normalized so norm with val only
                     enddo
                  enddo
               enddo!jsh
            enddo!ish
         enddo!jat
      enddo!iat

end

!! ------------------------------------------------------------------------
!  translate core shell string to integers
!! ------------------------------------------------------------------------

subroutine corestring(at,s)
      use cbascom
      implicit none
      integer at
      character*(*) s
      integer i,k,l,nsh

      if(index(s,'none').ne.0) return

      nsh=1
      k=0
      l=-1
      do i=1,len(s)/2
         l=l+2
         k=k+2
         if(s(k:k).eq.'s') cbas_lsh(nsh,at)=0
         if(s(k:k).eq.'p') cbas_lsh(nsh,at)=1
         if(s(k:k).eq.'d') cbas_lsh(nsh,at)=2
         if(s(l:l).eq.'1') cbas_npq(nsh,at)=1
         if(s(l:l).eq.'2') cbas_npq(nsh,at)=2
         if(s(l:l).eq.'3') cbas_npq(nsh,at)=3
         if(s(l:l).eq.'4') cbas_npq(nsh,at)=4
         if(s(l:l).eq.'5') cbas_npq(nsh,at)=5
         if(s(l:l).eq.'6') cbas_npq(nsh,at)=6
         nsh=nsh+1
      enddo

      nsh = nsh - 1
      cbas_nsh(at)=nsh

!     write(*,*) 'element',at
!     do i=1,nsh
!        write(*,*) cbas_npq(i,at),cbas_lsh(i,at)
!     enddo

end

!! ------------------------------------------------------------------------
!  setup core basis
!! ------------------------------------------------------------------------

subroutine setupcbas0(n,at)   
      use cbascom
      use  parcom
      implicit none
      integer n,at(n)

      integer i,il,iat,ish,iao,pr,maxpr,lladr(0:3),lladr2(0:3)
      data lladr  /1,3,6,10/
      data lladr2 /1,3,5, 7/
      character*20 cs(86)

      maxpr = 6

      cs='none'
      cbas_nsh=0
      cslexpo =0
      cslexpo =0
      clev    =0

      cs( 5:12)='1s'                       ! 2
      cs(13:30)='1s2s2p'                   !10
      cs(31:48)='1s2s2p3s3p3d'             !28
      cs(49:71)='1s2s2p3s3p3d4s4p4d'       !46
      cs(72:80)='1s2s2p3s3p3d4s4p4d'       !60-14=46 (f-core neglected)
      cs(81:86)='1s2s2p3s3p3d4s4p4d5s5p5d' !78-14=64    "       "

      do i=1,86
         call corestring(i,trim(cs(i))) ! set shells with n, l
      enddo

! Clementi's SZ Slater exponents from JCP 38 (1963), 2686 ibid 47 (1967), 1300.
      cslexpo(1,1) =  1.0000
      cslexpo(1,2) =  1.6250
      cslexpo(1,3) =  2.6906
      cslexpo(1,4) =  3.6848
      cslexpo(1,5) =  4.6795
      cslexpo(1,6) =  5.6727
      cslexpo(1,7) =  6.6651
      cslexpo(1,8) =  7.6579
      cslexpo(1,9) =  8.6501
      cslexpo(1,10)=  9.6421
      cslexpo(1,11)= 10.6259
      cslexpo(2,11)=  3.2857
      cslexpo(3,11)=  3.4409
      cslexpo(1,12)= 11.6089
      cslexpo(2,12)=  3.6960
      cslexpo(3,12)=  3.9129
      cslexpo(1,13)= 12.5910
      cslexpo(2,13)=  4.1068
      cslexpo(3,13)=  4.4817
      cslexpo(1,14)= 13.5745
      cslexpo(2,14)=  5.5100
      cslexpo(3,14)=  4.9725
      cslexpo(1,15)= 14.5578
      cslexpo(2,15)=  4.9125
      cslexpo(3,15)=  5.4806
      cslexpo(1,16)= 15.5409
      cslexpo(2,16)=  5.3144
      cslexpo(3,16)=  5.9855
      cslexpo(1,17)= 16.5239
      cslexpo(2,17)=  5.7152
      cslexpo(3,17)=  6.4966
      cslexpo(1,18)= 17.5075
      cslexpo(2,18)=  6.1152
      cslexpo(3,18)=  7.0041
      cslexpo(1,19)= 18.4895 
      cslexpo(2,19)=  6.5031
      cslexpo(3,19)=  7.5136 
      cslexpo(4,19)=  2.8933 
      cslexpo(5,19)=  2.5752 
      cslexpo(1,20)= 19.4730 
      cslexpo(2,20)=  6.8882
      cslexpo(3,20)=  8.0207 
      cslexpo(4,20)=  3.2005 
      cslexpo(5,20)=  2.8861 
      cslexpo(1,21)=  20.4566 
      cslexpo(2,21)=   7.2868 
      cslexpo(3,21)=   8.5273 
      cslexpo(1,22)=  21.4409 
      cslexpo(2,22)=   7.6883 
      cslexpo(3,22)=   9.0324 
      cslexpo(1,23)=  22.4256 
      cslexpo(2,23)=   8.0907 
      cslexpo(3,23)=   9.5364 
      cslexpo(1,24)=  23.4138 
      cslexpo(2,24)=   8.4919 
      cslexpo(3,24)=  10.0376 
      cslexpo(1,25)=  24.3957  
      cslexpo(2,25)=   8.8969 
      cslexpo(3,25)=  10.5420 
      cslexpo(1,26)=  25.3810 
      cslexpo(2,26)=   9.2995 
      cslexpo(3,26)=  11.0444 
      cslexpo(1,27)=  26.3668 
      cslexpo(2,27)=   9.7025 
      cslexpo(3,27)=  11.5462 
      cslexpo(1,28)=  27.3526 
      cslexpo(2,28)=  10.1063 
      cslexpo(3,28)=  12.0476 
      cslexpo(1,29)=  28.3386 
      cslexpo(2,29)=  10.5099 
      cslexpo(3,29)=  12.5485 
      cslexpo(1,30)=  29.3245 
      cslexpo(2,30)=  10.9140 
      cslexpo(3,30)=  13.0490 
      cslexpo(1,31)=  30.3094 
      cslexpo(2,31)=  11.2995 
      cslexpo(3,31)=  13.5454 
      cslexpo(4,31)=   5.6654 
      cslexpo(5,31)=   5.4012 
      cslexpo(6,31)=   5.0311 
      cslexpo(1,32)=  31.2937 
      cslexpo(2,32)=  11.6824 
      cslexpo(3,32)=  14.0411 
      cslexpo(4,32)=   5.9299 
      cslexpo(5,32)=   5.6712 
      cslexpo(6,32)=   5.4171 
      cslexpo(1,33)=  32.2783 
      cslexpo(2,33)=  12.0635 
      cslexpo(3,33)=  14.5368 
      cslexpo(4,33)=   6.1985 
      cslexpo(5,33)=   5.9499  
      cslexpo(6,33)=   5.7928 
      cslexpo(1,34)=  33.2622 
      cslexpo(2,34)=  12.4442 
      cslexpo(3,34)=  15.0326 
      cslexpo(4,34)=   6.4678 
      cslexpo(5,34)=   6.2350  
      cslexpo(6,34)=   6.1590 
      cslexpo(1,35)=  34.2471 
      cslexpo(2,35)=  12.8217 
      cslexpo(3,35)=  15.5282 
      cslexpo(4,35)=   6.7395 
      cslexpo(5,35)=   6.5236  
      cslexpo(6,35)=   6.5197 
      cslexpo(1,36)=  35.2316 
      cslexpo(2,36)=  13.1990 
      cslexpo(3,36)=  16.0235 
      cslexpo(4,36)=   7.0109 
      cslexpo(5,36)=   6.8114  
      cslexpo(6,36)=   6.8753 
      cslexpo(1,37)=36.2078 
      cslexpo(2,37)=13.5784 
      cslexpo(3,37)=16.5194 
      cslexpo(4,37)= 7.2809 
      cslexpo(5,37)= 7.1011
      cslexpo(6,37)= 7.2264 
      cslexpo(7,37)= 3.0970
      cslexpo(8,37)= 2.7202 
      cslexpo(1,38)=37.1911
      cslexpo(2,38)=13.9509
      cslexpo(3,38)=17.0152
      cslexpo(4,38)= 7.5546
      cslexpo(5,38)= 7.3892
      cslexpo(6,38)= 7.5754
      cslexpo(7,38)= 3.3611
      cslexpo(8,38)= 2.9830
      cslexpo(1,39)=38.1756 
      cslexpo(2,39)=14.3111
      cslexpo(3,39)=17.5016
      cslexpo(4,39)= 7.8505
      cslexpo(5,39)= 7.6975
      cslexpo(6,39)= 8.4657
      cslexpo(1,40)=39.1590 
      cslexpo(2,40)=14.6869
      cslexpo(3,40)=17.9964
      cslexpo(4,40)= 8.1205
      cslexpo(5,40)= 7.9485
      cslexpo(6,40)= 8.5223
      cslexpo(1,41)=40.1423 
      cslexpo(2,41)=15.0626
      cslexpo(3,41)=18.4911
      cslexpo(4,41)= 8.3905
      cslexpo(5,41)= 8.2184
      cslexpo(6,41)= 8.7847
      cslexpo(1,42)=41.1256 
      cslexpo(2,42)=15.4384
      cslexpo(3,42)=18.9859
      cslexpo(4,42)= 8.6605
      cslexpo(5,42)= 8.4912 
      cslexpo(6,42)= 9.0761
      cslexpo(1,43)=42.1090
      cslexpo(2,43)=15.8141
      cslexpo(3,43)=19.4704
      cslexpo(4,43)= 8.9304
      cslexpo(5,43)= 8.7947
      cslexpo(6,43)= 9.4510
      cslexpo(1,44)=43.0923
      cslexpo(2,44)=16.1899
      cslexpo(3,44)=19.9754
      cslexpo(4,44)= 9.2004
      cslexpo(5,44)= 9.0844
      cslexpo(6,44)= 9.7981
      cslexpo(1,45)=44.0756
      cslexpo(2,45)=16.5656
      cslexpo(3,45)=20.4702
      cslexpo(4,45)= 9.4704
      cslexpo(5,45)= 9.3724
      cslexpo(6,45)=10.1478
      cslexpo(1,46)=45.0589
      cslexpo(2,46)=16.9414
      cslexpo(3,46)=20.9650
      cslexpo(4,46)= 9.7404
      cslexpo(5,46)= 9.6616
      cslexpo(6,46)=10.4989
      cslexpo(1,47)=46.0423
      cslexpo(2,47)=17.3171
      cslexpo(3,47)=21.4597
      cslexpo(4,47)=10.0104
      cslexpo(5,47)= 9.9476
      cslexpo(6,47)=10.8503
      cslexpo(1,48)=47.0256
      cslexpo(2,48)=17.6929
      cslexpo(3,48)=21.9545
      cslexpo(4,48)=10.2804
      cslexpo(5,48)=10.2305
      cslexpo(6,48)=11.2023
      cslexpo(1,49)=48.0097
      cslexpo(2,49)=18.0618
      cslexpo(3,49)=22.4490
      cslexpo(4,49)=10.5436
      cslexpo(5,49)=10.5069
      cslexpo(6,49)=11.5594
      cslexpo(7,49)= 5.4403
      cslexpo(8,49)= 5.0922
      cslexpo(9,49)= 4.2354
      cslexpo(1,50)=48.9920
      cslexpo(2,50)=18.4297
      cslexpo(3,50)=22.9427
      cslexpo(4,50)=10.8066
      cslexpo(5,50)=10.7844
      cslexpo(6,50)=11.9139
      cslexpo(7,50)= 5.6645
      cslexpo(8,50)= 5.3163
      cslexpo(9,50)= 4.4925
      cslexpo(1,51)=49.9744
      cslexpo(2,51)=18.7977
      cslexpo(3,51)=23.4363
      cslexpo(4,51)=11.0697
      cslexpo(5,51)=11.0613
      cslexpo(6,51)=12.2666
      cslexpo(7,51)= 5.8859
      cslexpo(8,51)= 5.5453
      cslexpo(9,51)= 4.7436
      cslexpo(1,52)=50.9568
      cslexpo(2,52)=19.1656
      cslexpo(3,52)=23.9300
      cslexpo(4,52)=11.3327
      cslexpo(5,52)=11.3363
      cslexpo(6,52)=12.6131
      cslexpo(7,52)= 6.1021
      cslexpo(8,52)= 5.7805
      cslexpo(9,52)= 4.9900
      cslexpo(1,53)=51.9391 
      cslexpo(2,53)=19.5335 
      cslexpo(3,53)=24.4237 
      cslexpo(4,53)=11.5958
      cslexpo(5,53)=11.6138
      cslexpo(6,53)=12.9669
      cslexpo(7,53)= 6.3243
      cslexpo(8,53)= 6.0074 
      cslexpo(9,53)= 5.2335 
      cslexpo(1,54)=52.9215 
      cslexpo(2,54)=19.9015 
      cslexpo(3,54)=24.9173 
      cslexpo(4,54)=11.8588
      cslexpo(5,54)=11.8892
      cslexpo(6,54)=13.3156
      cslexpo(7,54)= 6.5432
      cslexpo(8,54)= 6.2393 
      cslexpo(9,54)= 5.4733 
      cslexpo(1,55)=53.9043
      cslexpo(2,55)=20.2558
      cslexpo(3,55)=25.4098
      cslexpo(4,55)=12.1258
      cslexpo(5,55)=12.1926
      cslexpo(6,55)=13.6602
      cslexpo(7,55)= 6.7606
      cslexpo(8,55)= 6.4644
      cslexpo(9,55)= 5.7096
      cslexpo(10,55)=3.0889
      cslexpo(11,55)=2.7302
      cslexpo(1, 56)=54.8861
      cslexpo(2, 56)=20.6234
      cslexpo(3, 56)=25.9048
      cslexpo(4, 56)=12.3852
      cslexpo(5, 56)=12.4388
      cslexpo(6, 56)=14.0081
      cslexpo(7, 56)= 6.9800
      cslexpo(8, 56)= 6.7008
      cslexpo(9, 56)= 5.9460
      cslexpo(10,56)= 3.3239
      cslexpo(11,56)= 2.9601
      cslexpo(1, 72)=70.6016
      cslexpo(2, 72)=26.5949
      cslexpo(3, 72)=33.7994
      cslexpo(4, 72)=16.7705
      cslexpo(5, 72)=16.9944
      cslexpo(6, 72)=19.4766
      cslexpo(7, 72)= 9.7443
      cslexpo(8, 72)= 9.4824
      cslexpo(9, 72)= 8.8810
      cslexpo(1, 73)=71.5837
      cslexpo(2, 73)=26.9649
      cslexpo(3, 73)=34.2932
      cslexpo(4, 73)=17.0305
      cslexpo(5, 73)=17.2668
      cslexpo(6, 73)=19.8137
      cslexpo(7, 73)= 9.9397
      cslexpo(8, 73)= 9.6837
      cslexpo(9, 73)= 9.0810
      cslexpo(1, 74)=72.5657
      cslexpo(2, 74)=27.3349
      cslexpo(3, 74)=34.7871
      cslexpo(4, 74)=17.2900
      cslexpo(5, 74)=17.5392
      cslexpo(6, 74)=20.1508
      cslexpo(7, 74)=10.1397
      cslexpo(8, 74)= 9.8871
      cslexpo(9, 74)= 9.2933
      cslexpo(1, 75)=73.5478
      cslexpo(2, 75)=27.7049
      cslexpo(3, 75)=35.2810
      cslexpo(4, 75)=17.5495
      cslexpo(5, 75)=17.8115
      cslexpo(6, 75)=20.4849
      cslexpo(7, 75)=10.3391
      cslexpo(8, 75)=10.0933
      cslexpo(9, 75)= 9.5136
      cslexpo(1, 76)=74.5299
      cslexpo(2, 76)=28.0749
      cslexpo(3, 76)=35.7749
      cslexpo(4, 76)=17.8091
      cslexpo(5, 76)=18.0839
      cslexpo(6, 76)=20.8249
      cslexpo(7, 76)=10.5238
      cslexpo(8, 76)=10.2860
      cslexpo(9, 76)= 9.7145
      cslexpo(1, 77)=75.5119
      cslexpo(2, 77)=28.4449
      cslexpo(3, 77)=36.2688
      cslexpo(4, 77)=18.0686
      cslexpo(5, 77)=18.3563
      cslexpo(6, 77)=21.1620
      cslexpo(7, 77)=10.7120
      cslexpo(8, 77)=10.4785
      cslexpo(9, 77)= 9.9343
      cslexpo(1, 78)=76.4940
      cslexpo(2, 78)=28.8149
      cslexpo(3, 78)=36.7627
      cslexpo(4, 78)=18.3281
      cslexpo(5, 78)=18.6287
      cslexpo(6, 78)=21.4991
      cslexpo(7, 78)=10.9097
      cslexpo(8, 78)=10.6826
      cslexpo(9, 78)=10.1575
      cslexpo(1, 79)=77.4761
      cslexpo(2, 79)=29.1849
      cslexpo(3, 79)=37.2566
      cslexpo(4, 79)=18.5876
      cslexpo(5, 79)=18.9010
      cslexpo(6, 79)=21.8361
      cslexpo(7, 79)=11.1033
      cslexpo(8, 79)=10.8867
      cslexpo(9, 79)=10.3820
      cslexpo(1, 80)=78.4581
      cslexpo(2, 80)=29.5547
      cslexpo(3, 80)=37.7505
      cslexpo(4, 80)=18.8471
      cslexpo(5, 80)=19.1734
      cslexpo(6, 80)=22.1732
      cslexpo(7, 80)=11.3112
      cslexpo(8, 80)=11.1015
      cslexpo(9, 80)=10.6170
      cslexpo(1, 81)=79.4409
      cslexpo(2, 81)=29.8421
      cslexpo(3, 81)=38.2431
      cslexpo(4, 81)=19.1397 
      cslexpo(5, 81)=19.4555
      cslexpo(6, 81)=22.5114
      cslexpo(7, 81)=11.5197 
      cslexpo(8, 81)=11.3042
      cslexpo(9, 81)=10.8472
      cslexpo(10,81)= 5.8244
      cslexpo(11,81)= 5.4177
      cslexpo(12,81)= 4.4050
      cslexpo(1, 82)=80.4195
      cslexpo(2, 82)=30.2150
      cslexpo(3, 82)=38.7383
      cslexpo(4, 82)=19.3841
      cslexpo(5, 82)=19.7165
      cslexpo(6, 82)=22.8489
      cslexpo(7, 82)=11.7232
      cslexpo(8, 82)=11.5084
      cslexpo(9, 82)=11.0799
      cslexpo(10,82)= 6.0263
      cslexpo(11,82)= 5.6060
      cslexpo(12,82)= 4.6304
      cslexpo(1, 83)=81.3982
      cslexpo(2, 83)=30.5880
      cslexpo(3, 83)=39.2335
      cslexpo(4, 83)=19.6285
      cslexpo(5, 83)=19.9774
      cslexpo(6, 83)=23.1805
      cslexpo(7, 83)=11.9268
      cslexpo(8, 83)=11.7126
      cslexpo(9, 83)=11.3098
      cslexpo(10,83)= 6.2058
      cslexpo(11,83)= 5.8042
      cslexpo(12,83)= 4.8488
      cslexpo(1, 84)=82.3768
      cslexpo(2, 84)=30.9609
      cslexpo(3, 84)=39.7286
      cslexpo(4, 84)=19.8729
      cslexpo(5, 84)=20.2383
      cslexpo(6, 84)=23.5240
      cslexpo(7, 84)=12.1304
      cslexpo(8, 84)=11.9168
      cslexpo(9, 84)=11.9168
      cslexpo(10,84)= 6.4046
      cslexpo(11,84)= 6.0049 
      cslexpo(12,84)= 5.0608  
      cslexpo(1, 85)=83.3554
      cslexpo(2, 85)=31.3338
      cslexpo(3, 85)=40.2238
      cslexpo(4, 85)=20.1173 
      cslexpo(5, 85)=20.4992
      cslexpo(6, 85)=23.8615
      cslexpo(7, 85)=12.3339
      cslexpo(8, 85)=12.1210
      cslexpo(9, 85)=11.7624
      cslexpo(10,85)= 6.5867 
      cslexpo(11,85)= 6.2080 
      cslexpo(12,85)= 5.2678
      cslexpo(1,86) =84.3341 
      cslexpo(2,86) =31.7068 
      cslexpo(3,86) =40.7190 
      cslexpo(4,86) =20.3617
      cslexpo(5,86) =20.7602
      cslexpo(6,86) =24.1991
      cslexpo(7,86) =12.5375
      cslexpo(8,86) =12.3253
      cslexpo(9,86) =11.9857
      cslexpo(10,86)= 6.7786
      cslexpo(11,86)= 6.3942
      cslexpo(12,86)= 5.4706

!     near HF levels
      clev(1, 1) =  -0.5000000
      clev(1, 2) =  -0.9179508
      clev(1, 3) =  -2.4777262
      clev(1, 4) =  -4.7863203
      clev(1, 5) =  -7.6940338
      clev(1, 6) = -11.3255170
      clev(1, 7) = -15.6290579
      clev(1, 8) = -20.6686540
      clev(1, 9) = -26.3827566
      clev(1,10) = -32.7724138
      clev(1,11) = -40.4782393
      clev(2,11) =  -2.7968751 
      clev(3,11) =  -1.5180133 
      clev(1,12) = -49.0801076
      clev(2,12) =  -3.8127947 
      clev(3,12) =  -2.3310958 
      clev(1,13) = -58.5009508
      clev(2,13) =  -4.9106092 
      clev(3,13) =  -3.2182410 
      clev(1,14) = -68.8123721
      clev(2,14) =  -6.1564677 
      clev(3,14) =  -4.2559826 
      clev(1,15) = -79.9696269
      clev(2,15) =  -7.5110229 
      clev(3,15) =  -5.4008827 
      clev(1,16) = -92.0043547
      clev(2,16) =  -9.0042114 
      clev(3,16) =  -6.6824236 
      clev(1,17) =-104.8843196
      clev(2,17) = -10.6073995
      clev(3,17) =  -8.0721369
      clev(1,18) =-118.6102421
      clev(2,18) = -12.3220674
      clev(3,18) =  -9.5713692
      clev(1,19) =-133.5237406
      clev(2,19) = -14.4775370
      clev(3,19) = -11.4944223
      clev(4,19) =  -1.7405194
      clev(5,19) =  -.94265914
      clev(1,20) =-149.3613693
      clev(2,20) = -16.8210000
      clev(3,20) = -13.6277555
      clev(4,20) =  -2.2444052
      clev(5,20) =  -1.3398417
      clev(1,21) =-165.8969001
      clev(2,21) = -19.0780965
      clev(3,21) = -15.6659817
      clev(1,22) =-183.2692614
      clev(2,22) = -21.4197345
      clev(3,22) = -17.7883022
      clev(1,23) =-201.4987403
      clev(2,23) = -23.8706785
      clev(3,23) = -20.0188512
      clev(1,24) =-220.5874573
      clev(2,24) = -26.4342375
      clev(3,24) = -22.3607530
      clev(1,25) =-240.5330490
      clev(2,25) = -29.1089634 
      clev(3,25) = -24.8121004
      clev(1,26) =-261.3667440
      clev(2,26) = -31.9280372
      clev(3,26) = -27.4067383
      clev(1,27) =-283.0575937
      clev(2,27) = -34.8591867
      clev(3,27) = -30.1116132
      clev(1,28) =-305.6097805
      clev(2,28) = -37.9068679
      clev(3,28) = -32.9314426
      clev(1,29) =-328.7838501
      clev(2,29) = -40.8100313
      clev(3,29) = -35.6095546
      clev(1,30) =-353.2923720
      clev(2,30) = -44.3467866
      clev(3,30) = -38.9107605
      clev(4,30) =  -5.6215957
      clev(5,30) =  -3.8244276
      clev(1,31) =-378.8176563
      clev(2,31) = -48.1674004
      clev(3,31) = -42.4930214 
      clev(4,31) =  -6.3935791 
      clev(5,31) =  -4.4813393 
      clev(6,31) =  -1.1921173 
      clev(1,32) =-405.2355329
      clev(2,32) = -52.1395939
      clev(3,32) = -46.2263666
      clev(4,32) =  -7.1797990 
      clev(5,32) =  -5.1513892 
      clev(6,32) =  -1.6202589 
      clev(1,33) =-432.5779002
      clev(2,33) = -56.2999996
      clev(3,33) = -50.1448668
      clev(4,33) =  -8.0195796 
      clev(5,33) =  -5.8715542 
      clev(6,33) =  -2.0988900 
      clev(1,34) =-460.8590130
      clev(2,34) = -60.6592102
      clev(3,34) = -54.2601596
      clev(4,34) =  -8.9223639 
      clev(5,34) =  -6.6526355 
      clev(6,34) =  -2.6359437 
      clev(1,35) =-490.0518044
      clev(2,35) = -65.1903812
      clev(3,35) = -58.5455164
      clev(4,35) =  -9.8623357 
      clev(5,35) =  -7.4694575 
      clev(6,35) =  -3.2064812 
      clev(1,36) =-520.1568748
      clev(2,36) = -69.8936596
      clev(3,36) = -63.0011463
      clev(4,36) = -10.8401409 
      clev(5,36) =  -8.3229134 
      clev(6,36) =  -3.8115953 
      clev(1,37) = -547.980732
      clev(2,37) =  -74.693343
      clev(3,37) =  -67.825593
      clev(4,37) =  -12.075449
      clev(5,37) =   -9.475439
      clev(6,37) =   -4.721253
      clev(7,37) =   -1.518441
      clev(8,37) =   -0.810443
      clev(1,38) = -579.760364 
      clev(2,38) =  -79.980107
      clev(3,38) =  -72.890869
      clev(4,38) =  -13.408268
      clev(5,38) =  -10.684558
      clev(6,38) =   -5.680669
      clev(7,38) =   -1.891043 
      clev(8,38) =   -1.099292
      clev(1,39)=  -612.239994
      clev(2,39)=   -85.216226
      clev(3,39)=   -77.906449
      clev(4,39)=   -14.551228
      clev(5,39)=   -11.704937
      clev(6,39)=    -6.455055
      clev(1,40)=  -645.621746 
      clev(2,40)=   -90.685052 
      clev(3,40)=   -83.168932 
      clev(4,40)=   -15.797094 
      clev(5,40)=   -12.828955 
      clev(6,40)=    -7.330992 
      clev(1,41)=  -679.938278 
      clev(2,41)=   -96.327380 
      clev(3,41)=   -88.599225 
      clev(4,41)=   -17.079740 
      clev(5,41)=   -13.987897 
      clev(6,41)=    -8.241586 
      clev(1,42)=   -715.14602
      clev(2,42)=   -102.15416
      clev(3,42)=   -94.198674
      clev(4,42)=   -18.401556
      clev(5,42)=   -15.181811
      clev(6,42)=    -9.186680
      clev(1,43)=  -751.332698
      clev(2,43)=  -108.238479
      clev(3,43)=  -100.067916
      clev(4,43)=   -19.862132
      clev(5,43)=   -16.523991
      clev(6,43)=   -10.282606
      clev(1,44)=  -788.313407
      clev(2,44)=  -114.416024
      clev(3,44)=  -106.026746
      clev(4,44)=   -21.277790
      clev(5,44)=   -17.803945
      clev(6,44)=   -11.321791
      clev(1,45)=  -826.232244 
      clev(2,45)=  -120.763025
      clev(3,45)=  -112.153807
      clev(4,45)=   -22.730477
      clev(5,45)=   -19.128054
      clev(6,45)=   -12.393816
      clev(1,46)=  -865.069606 
      clev(2,46)=  -127.281073 
      clev(3,46)=  -118.447098 
      clev(4,46)=   -24.220852  
      clev(5,46)=   -20.482548  
      clev(6,46)=   -13.508039  
      clev(1,47)=  -904.736227
      clev(2,47)=  -133.958897
      clev(3,47)=  -124.906322
      clev(4,47)=   -25.742792
      clev(5,47)=   -21.880871
      clev(6,47)=   -14.632828
      clev(1,48)=  -945.544783  
      clev(2,48)=  -141.029820  
      clev(3,48)=  -131.750662  
      clev(4,48)=   -27.521588   
      clev(5,48)=   -23.526435   
      clev(6,48)=   -16.024400   
      clev(7,48)=    -4.425086    
      clev(8,48)=    -3.050318    
      clev(1,49) =  -986.039127
      clev(2,49) =  -148.329602
      clev(3,49) =  -138.982025
      clev(4,49) =   -29.431590
      clev(5,49) =   -25.333321
      clev(6,49) =   -17.581263
      clev(7,49) =    -4.945096     
      clev(8,49) =    -3.504725    
      clev(9,49) =    -1.066632    
      clev(1,50) = -1028.462609
      clev(2,50) =  -155.822626
      clev(3,50) =  -146.261871
      clev(4,50) =   -31.394339
      clev(5,50) =   -27.165115
      clev(6,50) =   -19.156938
      clev(7,50) =    -5.482144     
      clev(8,50) =    -3.969944    
      clev(9,50) =    -1.377398
      clev(1,51) = -1071.885115
      clev(2,51) =  -163.525440
      clev(3,51) =  -153.721503
      clev(4,51) =   -33.415327
      clev(5,51) =   -29.048139
      clev(6,51) =   -20.782680
      clev(7,51) =    -6.026311      
      clev(8,51) =    -4.440078     
      clev(9,51) =    -1.692358
      clev(1,52) = -1116.072227 
      clev(2,52) =  -171.395384
      clev(3,52) =  -161.402426
      clev(4,52) =   -35.513193
      clev(5,52) =   -31.019088
      clev(6,52) =   -22.496886
      clev(7,52) =    -6.608889     
      clev(8,52) =    -4.949904    
      clev(9,52) =    -2.046883    
      clev(1,53) = -1177.186224
      clev(2,53) =  -180.949488
      clev(3,53) =  -169.660877
      clev(4,53) =   -37.934830
      clev(5,53) =   -33.122800
      clev(6,53) =   -24.286115
      clev(7,53) =    -7.244676    
      clev(8,53) =    -5.473820
      clev(9,53) =    -2.401629
      clev(1,54) = -1224.397652
      clev(2,54) =  -189.340333
      clev(3,54) =  -177.782921
      clev(4,54) =   -40.175973
      clev(5,54) =   -35.222101
      clev(6,54) =   -26.119243
      clev(7,54) =    -7.856564
      clev(8,54) =    -6.008750
      clev(9,54) =    -2.778258
      clev(1, 55)= -1262.576547  
      clev(2, 55)=  -196.743815
      clev(3, 55)=  -185.315843
      clev(4, 55)=   -42.418097
      clev(5, 55)=   -37.368061
      clev(6, 55)=   -28.085920
      clev(7, 55)=    -8.659781
      clev(8, 55)=    -6.739622
      clev(9, 55)=    -3.380502
      clev(10,55)=    -1.228517
      clev(11,55)=    -0.680752
      clev(1, 56)= -1311.680338
      clev(2, 56)=  -205.611145
      clev(3, 56)=  -193.816924
      clev(4, 56)=   -44.966581
      clev(5, 56)=   -39.744461
      clev(6, 56)=   -30.216149
      clev(7, 56)=    -9.514249     
      clev(8, 56)=    -7.506200
      clev(9, 56)=    -3.997419
      clev(10,56)=    -1.510172
      clev(11,56)=    -0.900567

      clev(1, 72)= -2213.645904
      clev(2, 72)=  -366.688200
      clev(3, 72)=  -348.985532
      clev(4, 72)=   -84.830032
      clev(5, 72)=   -76.585360
      clev(6, 72)=   -63.431921
      clev(7, 72)=   -17.630131
      clev(8, 72)=   -14.212758
      clev(9, 72)=    -8.753137
      clev(1, 73)= -2277.525256
      clev(2, 73)=  -378.519617
      clev(3, 73)=  -360.405013
      clev(4, 73)=   -87.994492
      clev(5, 73)=   -79.524638
      clev(6, 73)=   -66.160821
      clev(7, 73)=   -18.487062
      clev(8, 73)=   -14.956933
      clev(9, 73)=    -9.364571
      clev(1, 74)= -2342.330246
      clev(2, 74)=  -390.571409
      clev(3, 74)=  -372.032893
      clev(4, 74)=   -91.239797
      clev(5, 74)=   -82.538795
      clev(6, 74)=   -68.967486
      clev(7, 74)=   -19.373856
      clev(8, 74)=   -15.727283
      clev(9, 74)=   -10.002187
      clev(1, 75)= -2408.394783
      clev(2, 75)=  -402.936575
      clev(3, 75)=  -383.943938
      clev(4, 75)=   -94.654521
      clev(5, 75)=   -85.713161
      clev(6, 75)=   -71.938180
      clev(7, 75)=   -20.384014
      clev(8, 75)=   -16.618760
      clev(9, 75)=   -10.765282
      clev(1, 76)= -2475.611308
      clev(2, 76)=  -415.464765
      clev(3, 76)=  -395.986468
      clev(4, 76)=   -98.076899
      clev(5, 76)=   -88.882313
      clev(6, 76)=   -74.907027
      clev(7, 76)=   -21.344294
      clev(8, 76)=   -17.449423
      clev(9, 76)=   -11.469248
      clev(1, 77)= -2544.132094
      clev(2, 77)=  -428.228553
      clev(3, 77)=  -408.235971
      clev(4, 77)=  -101.578347
      clev(5, 77)=   -92.123900
      clev(6, 77)=   -77.949630
      clev(7, 77)=   -22.329022
      clev(8, 77)=   -18.304471
      clev(9, 77)=   -12.193095
      clev(1, 78)= -2613.096488
      clev(2, 78)=  -441.152576
      clev(3, 78)=  -420.681844
      clev(4, 78)=  -105.137875
      clev(5, 78)=   -95.429451
      clev(6, 78)=   -81.059442
      clev(7, 78)=   -23.332845
      clev(8, 78)=   -19.182828
      clev(9, 78)=   -12.943674
      clev(1, 79)= -2683.512776
      clev(2, 79)=  -454.317452
      clev(3, 79)=  -433.310346
      clev(4, 79)=  -108.768191
      clev(5, 79)=   -98.785511
      clev(6, 79)=   -84.218451
      clev(7, 79)=   -24.354862
      clev(8, 79)=   -20.062847
      clev(9, 79)=   -13.686616
      clev(1, 80)= -2754.748364
      clev(2, 80)=  -467.845242
      clev(3, 80)=  -446.308106
      clev(4, 80)=  -112.639779
      clev(5, 80)=  -102.381961
      clev(6, 80)=   -87.623964
      clev(7, 80)=   -25.578690
      clev(8, 80)=   -21.145629
      clev(9, 80)=   -14.641991
      clev(10,80)=    -4.193895
      clev(11,80)=    -2.751243
      clev(1, 81)= -2827.018225 
      clev(2, 81)=  -481.601770
      clev(3, 81)=  -459.453771
      clev(4, 81)=  -116.631912  
      clev(5, 81)=  -106.065182  
      clev(6, 81)=   -91.144789   
      clev(7, 81)=   -26.890410    
      clev(8, 81)=   -22.303972   
      clev(9, 81)=   -15.683246   
      clev(10,81)=    -4.639016    
      clev(11,81)=    -3.122645    
      clev(12,81)=    -0.985432
      clev(1, 82)= -2900.606343
      clev(2, 82)=  -495.578042
      clev(3, 82)=  -472.677922
      clev(4, 82)=  -120.730979  
      clev(5, 82)=  -109.792024   
      clev(6, 82)=   -94.769449   
      clev(7, 82)=   -28.248089    
      clev(8, 82)=   -23.484894   
      clev(9, 82)=   -16.768326    
      clev(10,82)=    -5.090038    
      clev(11,82)=    -3.492872  
      clev(12,82)=    -1.250144
      clev(1, 83)= -2975.353478
      clev(2, 83)=  -509.756597
      clev(3, 83)=  -486.071646
      clev(4, 83)=  -124.895174   
      clev(5, 83)=  -113.568438   
      clev(6, 83)=   -98.455627    
      clev(7, 83)=   -29.627671     
      clev(8, 83)=   -24.679241    
      clev(9, 83)=   -17.873706    
      clev(10,83)=    -5.539140     
      clev(11,83)=    -3.857371     
      clev(12,83)=    -1.512035     
      clev(1, 84)= -3050.371760
      clev(2, 84)=  -524.154675
      clev(3, 84)=  -499.687582
      clev(4, 84)=  -129.168353   
      clev(5, 84)=  -117.442762   
      clev(6, 84)=  -102.253755   
      clev(7, 84)=   -31.066816     
      clev(8, 84)=   -25.926638    
      clev(9, 84)=   -19.038207    
      clev(10,84)=    -6.019158     
      clev(11,84)=    -4.249033     
      clev(12,84)=    -1.804074     
      clev(1, 85)= -3126.802743 
      clev(2, 85)=  -538.771899
      clev(3, 85)=  -513.486907
      clev(4, 85)=  -133.511476
      clev(5, 85)=  -121.374707
      clev(6, 85)=  -106.120340
      clev(7, 85)=   -32.529169   
      clev(8, 85)=   -27.188944
      clev(9, 85)=   -20.223834
      clev(10,85)=    -6.499404
      clev(11,85)=    -4.637214
      clev(12,85)=    -2.096094
      clev(1, 86)= -3230.313904
      clev(2, 86)=  -556.914740
      clev(3, 86)=  -536.678708
      clev(4, 86)=  -138.423569 
      clev(5, 86)=  -128.673283 
      clev(6, 86)=  -110.703144 
      clev(7, 86)=   -33.922371  
      clev(8, 86)=   -29.492852 
      clev(9, 86)=   -21.333037 
      clev(10,86)=    -6.907222 
      clev(11,86)=    -5.226761 
      clev(12,86)=    -2.327901 

      allocate(corelist(n))
      ncorelist=0 
      do i=1,n
         if(abs(shell_cnf1(10,at(i))).gt.1d-6)then
            ncorelist=ncorelist+1
            corelist(ncorelist)=i
         endif
      enddo

      cncao=0
      cnsao=0
      cnpr =0
      do il=1,ncorelist
         i=corelist(il)
         iat = at(i)
         do ish=1,cbas_nsh(iat)
            cnsao = cnsao + lladr2(cbas_lsh(ish,iat))
            do iao=1,lladr(cbas_lsh(ish,iat))
               cncao = cncao + 1
               do pr=1,maxpr              
                  cnpr = cnpr + 1
               enddo
            enddo
         enddo
      enddo

      allocate(cprim_npr  (cncao))
      allocate(cprim_count(cncao))
      allocate(cprim_exp  (cnpr) )
      allocate(cprim_cnt  (cnpr) )
      allocate(caoshellc(13,n)   )
      allocate( aoshellc(13,n)   )  

end

!! ------------------------------------------------------------------------
!  setup core basis
!! ------------------------------------------------------------------------

subroutine setupcbas(n,at)
      use cbascom
      implicit none
      integer n,at(n)

      integer i,j,k,l,il,nn,iat,ish,iao,maxpr,lshell,pr
      integer lladr(0:3),lladr2(0:3)
      data lladr  /1,3,6,10/
      data lladr2 /1,3,5, 7/
      real*8 expo(6),cont(6)

      maxpr = 6

      k   = 0
      j   = 0
      cnpr= 0
      cncao=0
      cnsao=0
      do il=1,ncorelist
         i=corelist(il)
         iat = at(i)
         do ish=1,cbas_nsh(iat)
             caoshellc(ish,i)=k
              aoshellc(ish,i)=j
            lshell=cbas_lsh(ish,iat)
            k=k+lladr (lshell)          
            j=j+lladr2(lshell)               
            do iao=1,lladr2(lshell)              
               cnsao = cnsao + 1
            enddo
            call SETSTO6(cbas_npq(ish,iat),lshell+1,cslexpo(ish,iat),expo,cont)
            do iao=1,lladr(lshell)              
               cncao = cncao + 1
               cprim_npr(cncao)=maxpr            
               do pr=1,maxpr             
                  cnpr= cnpr+ 1
                  cprim_exp(cnpr) = expo(pr)
                  cprim_cnt(cnpr) = cont(pr)
                  if(lshell.eq.2.and.iao.gt.3) cprim_cnt(cnpr)=cprim_cnt(cnpr)*sqrt(3.0d0)
               enddo
            enddo
         enddo
      enddo

      k=0  
      do i=1,cncao
         cprim_count(i)=k
         k=k+cprim_npr(i)
      enddo   


      write(*,*) 'coreat ',ncorelist
      write(*,*) 'cncao  ',cncao
      write(*,*) 'cnsao  ',cnsao
      write(*,*) 'cnpr   ',cnpr 

end

!! ------------------------------------------------------------------------
!  setup core basis
!! ------------------------------------------------------------------------

subroutine modcbas(n,at,q)
      use cbascom
      implicit none
      integer n,at(n)
      real*8  q(n)

      integer i,j,k,l,nn,iat,ish,iao,npr,maxpr,lshell,pr
      integer lladr(0:3),lladr2(0:3)
      data lladr  /1,3,6,10/
      real*8 expo(6),cont(6),zeta

      maxpr = 6

      npr = 0
      do i=1,n   
         iat = at(i)
         do ish=1,cbas_nsh(iat)
            lshell=cbas_lsh(ish,iat)
            zeta  = cslexpo(ish,iat) ! *(1d0-q(i)*0.02d0)
            call SETSTO6(cbas_npq(ish,iat),lshell+1,zeta,expo,cont)
            do iao=1,lladr(lshell)              
               do pr=1,maxpr             
                  npr= npr+ 1
                  cprim_exp(npr) = expo(pr)
                  cprim_cnt(npr) = cont(pr)
                  if(lshell.eq.2.and.iao.gt.3) cprim_cnt(npr)=cprim_cnt(npr)*sqrt(3.0d0)
               enddo
            enddo
         enddo
      enddo

end

!! ------------------------------------------------------------------------
!  core core overlap for testing
!! ------------------------------------------------------------------------

subroutine csint2(nat,at,xyz,rab)
      use cbascom
      implicit none          
      integer, intent(in)  :: nat,at(nat)
      real*8, intent(in)   :: xyz(3,nat)
      real*8, intent(in)   :: rab(nat*(nat+1)/2)

      real*8 ss(6,6)
      real*8 tmp1,tmp2,tmp3,tmp4,intcut,s00,apb
      real*8 r1,r2,t1,t2,t3,t4,f,ci,cc,cj,alpi,rab2,ab,est
      real*8 rc(3),ra(3),rb(3),p(3),g(3),alpj,va(1),alpc,ccc
      integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kl,km,mi,mj,ij,lin,jshmax,ksh
      integer ip,jp,iat,jat,iatyp,jatyp,ish,jsh,icao,jcao,iao,jao,isao
      integer ishtyp,jshtyp,naoi,naoj,li,lj,iprim,jprim,iijj,iptyp,jptyp
      integer llao(0:3),llao2(0:3),itt(0:3),mli,mlj,abcnt
      data llao /1,3,6,10/
      data llao2/1,3,5,7 /
      data itt  /0,1,4,10/
      real*8 :: s(cnsao*(cnsao+1)/2)

      intcut=25. ! cut-off
      
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
            do ish=1,cbas_nsh(iatyp)
               ishtyp=cbas_lsh(ish,iatyp)
               icao=caoshellc(ish,iat)
               naoi=llao(ishtyp)
               iptyp=itt(ishtyp)
               jshmax=cbas_nsh(jatyp)
!              jshells
               do jsh=1,jshmax
                  jshtyp=cbas_lsh(jsh,jatyp)             
                  jcao=caoshellc(jsh,jat)
                  naoj=llao(jshtyp)
                  jptyp=itt(jshtyp)
                  ! we go through the primitives (beacause the screening is the same for all of them) 
                  ss=0.0d0
                  do ip=1,cprim_npr(icao+1)
                     iprim=ip+cprim_count(icao+1)
                     alpi=cprim_exp(iprim) ! exponent the same for each l component
                     do jp=1,cprim_npr(jcao+1)
                        jprim=jp+cprim_count(jcao+1)
                        alpj=cprim_exp(jprim) ! exponent the same for each l component
                        apb = alpi + alpj 
                        ab  =1.0d0/apb        
                        est =rab2*alpi*alpj*ab
                        if(est.gt.intcut) cycle

!                       now compute integrals for different components of i(e.g., px,py,pz)
! S and T
                        do mli=1,naoi
                           iprim=ip+cprim_count(icao+mli)
                           ci=cprim_cnt(iprim) ! coefficients NOT the same (contain CAO2SAO lin. comb. coefficients)
                           do mlj=1,naoj
                             jprim=jp+cprim_count(jcao+mlj)
                             cc=cprim_cnt(jprim)*ci
                             call propa_s  (nat,xyz,ra,rb,alpi,alpj,iptyp+mli,jptyp+mlj,va)
                             ss(mlj,mli)=ss(mlj,mli)+va(1)*cc
                           enddo ! mlj
                        enddo ! mli
                     enddo ! jp
                  enddo ! ip 
                  !transform from CAO to SAO
                  call dtrf2(ss,ishtyp,jshtyp)
                  do ii=1,llao2(ishtyp)
                     iao=ii+aoshellc(ish,iat)
                     do jj=1,llao2(jshtyp)
                        jao=jj+aoshellc(jsh,jat)
                        ij = lin(jao,iao)
                        s(ij)=ss(jj,ii)
                     enddo
                  enddo
               enddo!jsh
            enddo!ish
         enddo!jat
      enddo!iat

      ij = 0
      do i=1,cnsao
         ij = ij + i
         write(*,*) 1./sqrt(s(ij))
      enddo

end

