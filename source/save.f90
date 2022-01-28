subroutine int3gen12(nat,nao,at,xyz,rab,norm,v)
      use bascom
      use parcom
      implicit none          
      integer, intent(in)  :: nao,nat
      integer, intent(in)  :: at(nat)     
      real*8, intent(in)   :: xyz(3,nat)
      real*8, intent(in)   :: rab(nat*(nat+1)/2)
      real*8, intent(in)   :: norm(nao)
      real*8, intent(inout):: v(2*naux+1,nao*(nao+1)/2)

      real*8 cc,alpi,rab2,ab,est,alpc,ccc,ff
      real*8 ra(3),rb(3),rc(3),alpj,Xab,apb,pa,pb
      real*8 g(3),p(3),s00,rpq2,eabc,abc,s000,expe
      real*8 pindex(0:13)
      real*8 cccc(naux)   
      integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kl,km,mi,mj,ij,lin,jshmax
      integer ip,jp,iat,jat,iatyp,jatyp,ish,jsh,icao,jcao,iao,jao,isao,ksh
      integer ishtyp,jshtyp,naoi,naoj,li,lj,iprim,jprim,kadr,icase
      integer llao(0:3),llao2(0:3),mli,mlj,iaux,prop
      data llao /1,3,6,10/
      data llao2/1,3,5,7 /
      real*8 vv(6,6,2*naux+1), va2(6,6),va3(6,6)
      real*8, parameter :: SQRT_SRF = sqrt(139.94734662099890277426d0)
      real*8, parameter ::       PI = 3.1415926535897932384626433832795029d0
      real*8, parameter ::       PI2= 2.0d0/PI
      real*8, parameter :: aux_norm = SQRT_SRF/(2.0d0*sqrt(2.d0))               
      real*8, parameter ::      thr = 1.d-6 ! aggressive cut-off

      kadr=2*naux
      v=0
      do iat=1,nat
         ra(1:3)=xyz(1:3,iat)
         iatyp=at(iat)
         do jat=1,iat
            rb(1:3)=xyz(1:3,jat)
            jatyp=at(jat)
            rab2=rab(lin(iat,jat))**2               
            if(rab2.gt.200) cycle
            do ish=1,bas_nsh(iatyp)
               ishtyp=bas_lsh(ish,iatyp)
               icao=caoshell(ish,iat)
               naoi=llao(ishtyp)
               jshmax=bas_nsh(jatyp)
!              jshells
               do jsh=1,jshmax
                  jshtyp=bas_lsh(jsh,jatyp)             
                  jcao=caoshell(jsh,jat)
                  naoj=llao(jshtyp)
                  select case (ishtyp)
                    case (0)
                    select case (jshtyp)
                      case (0)
                      icase=1
                      case (1)
                      icase=2
                      case (2)
                      icase=3
                    end select
                    case (1)
                    select case (jshtyp)
                      case (0)
                      icase=4
                      case (1)
                      icase=5
                      case (2)
                      icase=6
                    end select
                    case (2)
                    select case (jshtyp)
                      case (0)
                      icase=7
                      case (1)
                      icase=8
                      case (2)
                      icase=9
                    end select
                  end select
                  vv=0
                  do ip=1,prim_npr(icao+1)
                     iprim=ip+prim_count(icao+1)
                     alpi=prim_exp(iprim) ! exponent 
                     do jp=1,prim_npr(jcao+1)
                        jprim=jp+prim_count(jcao+1)
                        alpj=prim_exp(jprim) ! exponent 
                        apb =alpi+alpj
                        ab  =1.0d0/apb
                        est =rab2*alpi*alpj*ab
                        cc  =prim_cnt(jprim)*prim_cnt(iprim)  ! contraction coeff
                        expe=exp(-est)                 
                        Xab =SQRT_SRF*cc*expe                 
                        if(abs(Xab).lt.thr) cycle     ! cut prim pair
!!!!!!!!!!!
! <i|s/r|j>
!!!!!!!!!!!
                        pa  =alpi/apb
                        pb  =alpj/apb
                        pindex(0) = 0.5d0*ab
                        pindex(1) = apb
                        pindex(2) = Xab
                        pindex(3) = pa*ra(1)+pb*rb(1)
                        pindex(4) = pa*ra(2)+pb*rb(2)
                        pindex(5) = pa*ra(3)+pb*rb(3)
                        pindex(6) = pindex(3)-ra(1)    
                        pindex(7) = pindex(4)-ra(2)    
                        pindex(8) = pindex(5)-ra(3)    
                        pindex(9) = pindex(3)-rb(1)    
                        pindex(10)= pindex(4)-rb(2)    
                        pindex(11)= pindex(5)-rb(3)    
                        pindex(12)= 2.*alpi            
                        pindex(13)= 2.*alpj            
                        s00=(pi*ab)**1.50d0*expe
                        p(:)=(alpi*ra(:)+alpj*rb(:))*ab
                        select case (icase)
                        case (1)
                             do ksh=1, naux        
                             rc(:)=xyz(:,aux_at(ksh))
                             alpc =aux_exp(ksh)   
                             ccc=cccc(ksh)
                             eabc= 1.0d0/(apb + alpc)
                             g(:)=(apb*p(:)+alpc*rc(:))*eabc
                             rpq2=(p(1)-rc(1))**2+(p(2)-rc(2))**2+(p(3)-rc(3))**2
                             s000= s00*exp(-(apb*alpc*eabc)*rpq2)*(apb*eabc)**1.50d0
                             call sss2e(pindex,alpc,ccc,rc,rpq2,va2)
                             call ss1e (ra,rb,g,eabc,s000,va3)
                             vv(:,:,ksh)     =vv(:,:,ksh)     +va2(:,:)
                             vv(:,:,ksh+naux)=vv(:,:,ksh+naux)+va3(:,:)*ccc*cc/aux_norm
                             enddo
                        case (2)
                             do ksh=1, naux        
                             rc(:)=xyz(:,aux_at(ksh))
                             alpc =aux_exp(ksh)   
                             ccc=cccc(ksh) 
                             eabc= 1.0d0/(apb + alpc)
                             g(:)=(apb*p(:)+alpc*rc(:))*eabc
                             rpq2=(p(1)-rc(1))**2+(p(2)-rc(2))**2+(p(3)-rc(3))**2
                             s000= s00*exp(-(apb*alpc*eabc)*rpq2)*(apb*eabc)**1.50d0
                             call sps2e(pindex,alpc,ccc,rc,rpq2,va2)
                             call sp1e (ra,rb,g,eabc,s000,va3)
                             vv(:,:,ksh)     =vv(:,:,ksh)     +va2(:,:)
                             vv(:,:,ksh+naux)=vv(:,:,ksh+naux)+va3(:,:)*ccc*cc/aux_norm
                             enddo
                        case (3)
                             do ksh=1, naux        
                             rc(:)=xyz(:,aux_at(ksh))
                             alpc =aux_exp(ksh)   
                             ccc=cccc(ksh)
                             eabc= 1.0d0/(apb + alpc)
                             g(:)=(apb*p(:)+alpc*rc(:))*eabc
                             rpq2=(p(1)-rc(1))**2+(p(2)-rc(2))**2+(p(3)-rc(3))**2
                             s000= s00*exp(-(apb*alpc*eabc)*rpq2)*(apb*eabc)**1.50d0
                             call sds2e(pindex,alpc,ccc,rc,rpq2,va2)
                             call sd1e (ra,rb,g,eabc,s000,va3)
                             vv(:,:,ksh)     =vv(:,:,ksh)     +va2(:,:)
                             vv(:,:,ksh+naux)=vv(:,:,ksh+naux)+va3(:,:)*ccc*cc/aux_norm
                             enddo
                        case (4)
                             do ksh=1, naux        
                             rc(:)=xyz(:,aux_at(ksh))
                             alpc =aux_exp(ksh)   
                             ccc=cccc(ksh)
                             eabc= 1.0d0/(apb + alpc)
                             g(:)=(apb*p(:)+alpc*rc(:))*eabc
                             rpq2=(p(1)-rc(1))**2+(p(2)-rc(2))**2+(p(3)-rc(3))**2
                             s000= s00*exp(-(apb*alpc*eabc)*rpq2)*(apb*eabc)**1.50d0
                             call pss2e(pindex,alpc,ccc,rc,rpq2,va2)
                             call ps1e (ra,rb,g,eabc,s000,va3)
                             vv(:,:,ksh)     =vv(:,:,ksh)     +va2(:,:)
                             vv(:,:,ksh+naux)=vv(:,:,ksh+naux)+va3(:,:)*ccc*cc/aux_norm
                             enddo
                        case (5)
                             do ksh=1, naux        
                             rc(:)=xyz(:,aux_at(ksh))
                             alpc =aux_exp(ksh)   
                             ccc=cccc(ksh)
                             eabc= 1.0d0/(apb + alpc)
                             g(:)=(apb*p(:)+alpc*rc(:))*eabc
                             rpq2=(p(1)-rc(1))**2+(p(2)-rc(2))**2+(p(3)-rc(3))**2
                             s000= s00*exp(-(apb*alpc*eabc)*rpq2)*(apb*eabc)**1.50d0
                             call pps2e(pindex,alpc,ccc,rc,rpq2,va2)
                             call pp1e (ra,rb,g,eabc,s000,va3)
                             vv(:,:,ksh)     =vv(:,:,ksh)     +va2(:,:)
                             vv(:,:,ksh+naux)=vv(:,:,ksh+naux)+va3(:,:)*ccc*cc/aux_norm
                             enddo
                        case (6)
                             do ksh=1, naux        
                             rc(:)=xyz(:,aux_at(ksh))
                             alpc =aux_exp(ksh)   
                             ccc=cccc(ksh)
                             eabc= 1.0d0/(apb + alpc)
                             g(:)=(apb*p(:)+alpc*rc(:))*eabc
                             rpq2=(p(1)-rc(1))**2+(p(2)-rc(2))**2+(p(3)-rc(3))**2
                             s000= s00*exp(-(apb*alpc*eabc)*rpq2)*(apb*eabc)**1.50d0
                             call pds2e(pindex,alpc,ccc,rc,rpq2,va2)
                             call pd1e (ra,rb,g,eabc,s000,va3)
                             vv(:,:,ksh)     =vv(:,:,ksh)     +va2(:,:)
                             vv(:,:,ksh+naux)=vv(:,:,ksh+naux)+va3(:,:)*ccc*cc/aux_norm
                             enddo
                        case (7)
                             do ksh=1, naux        
                             rc(:)=xyz(:,aux_at(ksh))
                             alpc =aux_exp(ksh)   
                             ccc=cccc(ksh)
                             eabc= 1.0d0/(apb + alpc)
                             g(:)=(apb*p(:)+alpc*rc(:))*eabc
                             rpq2=(p(1)-rc(1))**2+(p(2)-rc(2))**2+(p(3)-rc(3))**2
                             s000= s00*exp(-(apb*alpc*eabc)*rpq2)*(apb*eabc)**1.50d0
                             call dss2e(pindex,alpc,ccc,rc,rpq2,va2)
                             call ds1e (ra,rb,g,eabc,s000,va3)
                             vv(:,:,ksh)     =vv(:,:,ksh)     +va2(:,:)
                             vv(:,:,ksh+naux)=vv(:,:,ksh+naux)+va3(:,:)*ccc*cc/aux_norm
                             enddo
                        case (8)
                             do ksh=1, naux        
                             rc(:)=xyz(:,aux_at(ksh))
                             alpc =aux_exp(ksh)   
                             ccc=cccc(ksh)
                             eabc= 1.0d0/(apb + alpc)
                             g(:)=(apb*p(:)+alpc*rc(:))*eabc
                             rpq2=(p(1)-rc(1))**2+(p(2)-rc(2))**2+(p(3)-rc(3))**2
                             s000= s00*exp(-(apb*alpc*eabc)*rpq2)*(apb*eabc)**1.50d0
                             call dps2e(pindex,alpc,ccc,rc,rpq2,va2)
                             call dp1e (ra,rb,g,eabc,s000,va3)
                             vv(:,:,ksh)     =vv(:,:,ksh)     +va2(:,:)
                             vv(:,:,ksh+naux)=vv(:,:,ksh+naux)+va3(:,:)*ccc*cc/aux_norm
                             enddo
                        case (9)
                             do ksh=1, naux        
                             rc(:)=xyz(:,aux_at(ksh))
                             alpc =aux_exp(ksh)   
                             ccc=cccc(ksh)
                             eabc= 1.0d0/(apb + alpc)
                             g(:)=(apb*p(:)+alpc*rc(:))*eabc
                             rpq2=(p(1)-rc(1))**2+(p(2)-rc(2))**2+(p(3)-rc(3))**2
                             s000= s00*exp(-(apb*alpc*eabc)*rpq2)*(apb*eabc)**1.50d0
                             call dds2e(pindex,alpc,ccc,rc,rpq2,va2)
                             call dd1e (ra,rb,g,eabc,s000,va3)
                             vv(:,:,ksh)     =vv(:,:,ksh)     +va2(:,:)
                             vv(:,:,ksh+naux)=vv(:,:,ksh+naux)+va3(:,:)*ccc*cc/aux_norm
                             enddo
                        end select

                     enddo ! jp
                  enddo ! ip 

                  !transform from CAO to SAO
                  if(ishtyp.eq.2.or.jshtyp.eq.2)then
                    do prop=1,kadr 
                     call dtrf3(vv(:,:,prop),ishtyp,jshtyp)
                    enddo
                  endif
                  do ii=1,llao2(ishtyp)
                     iao=ii+aoshell(ish,iat)
                     do jj=1,llao2(jshtyp)
                        jao=jj+aoshell(jsh,jat)
                        ij = lin(iao,jao)   
                        ff=norm(iao)*norm(jao)
                        do prop=1,kadr 
                           v(prop,ij)=vv(jj,ii,prop)*ff
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

end subroutine int3gen12


