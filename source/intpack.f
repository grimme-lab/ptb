c=======================================================================
c cartesian gaussian functions (6d,10f...)
c iff :
c s,px, py pz, dx**2 dy**2 dz**2 dxy dxz dyz
c 1 2   3   4   5     6     7     8   9  10
c fxxx, fyyy, fzzz, fxxy, fxxz, fyyx, fyyz, fxzz, fyzz, fxyz
c   11   12    13    14    15    16    17    18   19    20
c a, are aufpunkte
c nt : # of returns
c va : integral
c=======================================================================

      subroutine propa_stv(nuc,xyz,a,b,etaij4,etakl4,iff1,iff2,va)   
      implicit real*8(a-h,o-z)                                                  
c aufpunkte,ref point,intarray
      real*8 a(3),b(3),va(nuc+2),xyz(3,nuc)                                     
c local
      common /abfunc/ ra(3),rb(3),ga,gb,ia,ib
      dimension d(3),dd(84),v(nuc+2)
      dimension e(3),aa(20),bb(20)
      dimension lin(84),min(84),nin(84)                                         
      data lin/0,1,0,0,2,0,0,1,1,0,3,0,0,2,2,1,0,1,0,1,4,0,0,3,3,1,0,1,0        
     1 ,2,2,0,2,1,1,5,0,0,3,3,2,2,0,0,4,4,1,0,0,1,1,3,1,2,2,1,6,0,0,3,          
     2 3,0,5,5,1,0,0,1,4,4,2,0,2,0,3,3,1,2,2,1,4,1,1,2/,min/0,0,1,0,0,          
     3 2,0,1,0,1,0,3,0,1,0,2,2,0,1,1,0,4,0,1,0,3,3,0,1,2,0,2,1,2,1,0,5,         
     4 0,2,0,3,0,3,2,1,0,4,4,1,0,1,1,3,2,1,2,0,6,0,3,0,3,1,0,0,1,5,5,2,0        
     5,0,2,4,4,2,1,3,1,3,2,1,4,1,2/,nin/0,0,0,1,0,0,2,0,1,1,0,0,3,0,1,0,        
     6 1,2,2,1,0,0,4,0,1,0,1,3,3,0,2,2,1,1,2,0,0,5,0,2,0,3,2,3,0,1,0,1,         
     7 4,4,3,1,1,1,2,2,0,0,6,0,3,3,0,1,5,5,1,0,0,2,4,4,0,2,1,2,2,3,1,3,         
     8 1,1,4,2/                                                                 

c --- a,b are centres of gaussians                                              
      ra = a
      rb = b
c --- ga,gb are their exponents                                                 
      ga=etaij4                                                              
      gb=etakl4                                                              
      aa = 0
      bb = 0
      aa(iff1)=1.0d0
      bb(iff2)=1.0d0
      ia=iff1                                                                   
      ib=iff2                                                                   
c --- ia,ib are the canonical indices for monom                                 
c --- apply product theorem                                                     
      cij=0.0d0
      ckl=0.0d0
      call divpt(a,etaij4,b,etakl4,cij,ckl,e,gama,efact)                  

c --- calculate cartesian prefactor for first gaussian                          
      call rhftce(aa,a,e,iff1)                                                  
c --- calculate cartesian prefactor for second gaussian                         
      call rhftce(bb,b,e,iff2)                                                  
c --- form their product                                                        
      call prod(aa,bb,dd,iff1,iff2)                                             
      va = 0
c ----- e is center of product gaussian with exponent gama                      
c          aname represents an external function                                
      if(iff1.gt.10.or.iff2.gt.10) goto 110                                     
      if(iff1.gt. 4.or.iff2.gt. 4) goto 120                                       
      if(iff1.gt. 1.or.iff2.gt. 1) goto 130                                       
c           s - s                                                               
      d = e 
      call opad1(lin(1),min(1),nin(1),gama,v(1),d)                                 
      call opap4(lin(1),min(1),nin(1),gama,v(2),d)                                 
!$OMP PARALLEL PRIVATE ( j, d ) SHARED ( xyz,e,v,gama,lin,min,nin )
!$OMP DO
      do j=1,nuc
         d(1:3) = e(1:3) - xyz(1:3,j)
         call opaa0(lin(1),min(1),nin(1),gama,v(j+2),d)                                 
      enddo
!$OMP END DO
!$OMP END PARALLEL      
      do j=1,nuc+2                                                          
         va(j)=dd(1)*v(j)*efact                                                   
      enddo
      return                                                                    
  130 continue                                                                  
      if(iff1.ne.1.and.iff2.ne.1) goto 131                                      
c           s - p                                                               
      do i=1,4                                                               
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         d = e 
         call opad1(lin(i),min(i),nin(i),gama,v(1),d)                                 
         call opap4(lin(i),min(i),nin(i),gama,v(2),d)                                 
!$OMP PARALLEL PRIVATE ( j, d ) SHARED ( xyz,e,v,gama,lin,min,nin )
!$OMP DO
         do j=1,nuc
            d(1:3) = e(1:3) - xyz(1:3,j)
            call opaa0(lin(i),min(i),nin(i),gama,v(j+2),d)                                 
         enddo
!$OMP END DO
!$OMP END PARALLEL      
         do j=1,nuc+2                                                          
            va(j)=va(j)+dd(i)*v(j)*efact                                                   
         enddo
      enddo
      return                                                                    
  131 continue                                                                  
c           p - p                                                               
      do i=1,10                                                              
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         d = e 
         call opad1(lin(i),min(i),nin(i),gama,v(1),d)                                 
         call opap4(lin(i),min(i),nin(i),gama,v(2),d)                                 
!$OMP PARALLEL PRIVATE ( j, d ) SHARED ( xyz,e,v,gama,lin,min,nin )
!$OMP DO
         do j=1,nuc
            d(1:3) = e(1:3) - xyz(1:3,j)
            call opaa0(lin(i),min(i),nin(i),gama,v(j+2),d)                                 
         enddo
!$OMP END DO
!$OMP END PARALLEL      
         do j=1,nuc+2                                                          
            va(j)=va(j)+dd(i)*v(j)*efact                                                   
         enddo
      enddo
      return                                                                    
  120 continue                                                                  
      if(iff1.ne.1.and.iff2.ne.1) goto 121                                      
c           s - d                                                               
      do i=1,10                                                              
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         d = e 
         call opad1(lin(i),min(i),nin(i),gama,v(1),d)                                 
         call opap4(lin(i),min(i),nin(i),gama,v(2),d)                                 
!$OMP PARALLEL PRIVATE ( j, d ) SHARED ( xyz,e,v,gama,lin,min,nin )
!$OMP DO
         do j=1,nuc
            d(1:3) = e(1:3) - xyz(1:3,j)
            call opaa0(lin(i),min(i),nin(i),gama,v(j+2),d)                                 
         enddo
!$OMP END DO
!$OMP END PARALLEL      
         do j=1,nuc+2                                                          
            va(j)=va(j)+dd(i)*v(j)*efact                                                   
         enddo
      enddo
      return                                                                    
  121 continue                                                                  
      if(iff1.gt.4.and.iff2.gt.4) goto 122                                      
c           p - d                                                               
      do i=1,20                                                              
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         d = e 
         call opad1(lin(i),min(i),nin(i),gama,v(1),d)                                 
         call opap4(lin(i),min(i),nin(i),gama,v(2),d)                                 
!$OMP PARALLEL PRIVATE ( j, d ) SHARED ( xyz,e,v,gama,lin,min,nin )
!$OMP DO
         do j=1,nuc
            d(1:3) = e(1:3) - xyz(1:3,j)
            call opaa0(lin(i),min(i),nin(i),gama,v(j+2),d)                                 
         enddo
!$OMP END DO     
!$OMP END PARALLEL      
         do j=1,nuc+2                                                          
            va(j)=va(j)+dd(i)*v(j)*efact                                                   
         enddo
      enddo
      return                                                                    
  122 continue                                                                  
c           d - d                                                               
      do i=1,35                                                              
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         d = e 
         call opad1(lin(i),min(i),nin(i),gama,v(1),d)                                 
         call opap4(lin(i),min(i),nin(i),gama,v(2),d)                                 
!$OMP PARALLEL PRIVATE ( j, d ) SHARED ( xyz,e,v,gama,lin,min,nin )
!$OMP DO
         do j=1,nuc
            d(1:3) = e(1:3) - xyz(1:3,j)
            call opaa0(lin(i),min(i),nin(i),gama,v(j+2),d)                                 
         enddo
!$OMP END DO     
!$OMP END PARALLEL      
         do j=1,nuc+2                                                          
            va(j)=va(j)+dd(i)*v(j)*efact                                                   
         enddo
      enddo
      return                                                                    
  110 continue                                                                  
      stop 'no f-fkt.'

      end                                                                       

c=======================================================================
c cartesian gaussian functions (6d,10f...)
c iff :
c s,px, py pz, dx**2 dy**2 dz**2 dxy dxz dyz
c 1 2   3   4   5     6     7     8   9  10
c fxxx, fyyy, fzzz, fxxy, fxxz, fyyx, fyyz, fxzz, fyzz, fxyz
c   11   12    13    14    15    16    17    18   19    20
c a, are aufpunkte
c nt : # of returns
c va : integral
c=======================================================================

      subroutine propa_sv(nuc,xyz,a,b,etaij4,etakl4,iff1,iff2,va)   
      implicit real*8(a-h,o-z)                                                  
c aufpunkte,ref point,intarray
      real*8 a(3),b(3),va(nuc+1),xyz(3,nuc)                                     
c local
      common /abfunc/ ra(3),rb(3),ga,gb,ia,ib
      dimension d(3),dd(84),v(nuc+1)
      dimension e(3),aa(20),bb(20)
      dimension lin(84),min(84),nin(84)                                         
      data lin/0,1,0,0,2,0,0,1,1,0,3,0,0,2,2,1,0,1,0,1,4,0,0,3,3,1,0,1,0        
     1 ,2,2,0,2,1,1,5,0,0,3,3,2,2,0,0,4,4,1,0,0,1,1,3,1,2,2,1,6,0,0,3,          
     2 3,0,5,5,1,0,0,1,4,4,2,0,2,0,3,3,1,2,2,1,4,1,1,2/,min/0,0,1,0,0,          
     3 2,0,1,0,1,0,3,0,1,0,2,2,0,1,1,0,4,0,1,0,3,3,0,1,2,0,2,1,2,1,0,5,         
     4 0,2,0,3,0,3,2,1,0,4,4,1,0,1,1,3,2,1,2,0,6,0,3,0,3,1,0,0,1,5,5,2,0        
     5,0,2,4,4,2,1,3,1,3,2,1,4,1,2/,nin/0,0,0,1,0,0,2,0,1,1,0,0,3,0,1,0,        
     6 1,2,2,1,0,0,4,0,1,0,1,3,3,0,2,2,1,1,2,0,0,5,0,2,0,3,2,3,0,1,0,1,         
     7 4,4,3,1,1,1,2,2,0,0,6,0,3,3,0,1,5,5,1,0,0,2,4,4,0,2,1,2,2,3,1,3,         
     8 1,1,4,2/                                                                 

c --- a,b are centres of gaussians                                              
      ra = a
      rb = b
c --- ga,gb are their exponents                                                 
      ga=etaij4                                                              
      gb=etakl4                                                              
      aa = 0
      bb = 0
      aa(iff1)=1.0d0
      bb(iff2)=1.0d0
      ia=iff1                                                                   
      ib=iff2                                                                   
c --- ia,ib are the canonical indices for monom                                 
c --- apply product theorem                                                     
      cij=0.0d0
      ckl=0.0d0
      call divpt(a,etaij4,b,etakl4,cij,ckl,e,gama,efact)                  

c --- calculate cartesian prefactor for first gaussian                          
      call rhftce(aa,a,e,iff1)                                                  
c --- calculate cartesian prefactor for second gaussian                         
      call rhftce(bb,b,e,iff2)                                                  
c --- form their product                                                        
      call prod(aa,bb,dd,iff1,iff2)                                             
      va = 0
c ----- e is center of product gaussian with exponent gama                      
c          aname represents an external function                                
      if(iff1.gt.10.or.iff2.gt.10) goto 110                                     
      if(iff1.gt. 4.or.iff2.gt. 4) goto 120                                       
      if(iff1.gt. 1.or.iff2.gt. 1) goto 130                                       
c           s - s                                                               
      d = e 
      call opad1(lin(1),min(1),nin(1),gama,v(1),d)                                 
!$OMP PARALLEL PRIVATE ( j, d ) SHARED ( xyz,e,v,gama,lin,min,nin )
!$OMP DO
      do j=1,nuc
         d(1:3) = e(1:3) - xyz(1:3,j)
         call opaa0(lin(1),min(1),nin(1),gama,v(j+1),d)                                 
      enddo
!$OMP END DO
!$OMP END PARALLEL      
      do j=1,nuc+1                                                          
         va(j)=dd(1)*v(j)*efact                                                   
      enddo
      return                                                                    
  130 continue                                                                  
      if(iff1.ne.1.and.iff2.ne.1) goto 131                                      
c           s - p                                                               
      do i=1,4                                                               
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         d = e 
         call opad1(lin(i),min(i),nin(i),gama,v(1),d)                                 
!$OMP PARALLEL PRIVATE ( j, d ) SHARED ( xyz,e,v,gama,lin,min,nin )
!$OMP DO
         do j=1,nuc
            d(1:3) = e(1:3) - xyz(1:3,j)
            call opaa0(lin(i),min(i),nin(i),gama,v(j+1),d)                                 
         enddo
!$OMP END DO
!$OMP END PARALLEL      
         do j=1,nuc+1                                                          
            va(j)=va(j)+dd(i)*v(j)*efact                                                   
         enddo
      enddo
      return                                                                    
  131 continue                                                                  
c           p - p                                                               
      do i=1,10                                                              
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         d = e 
         call opad1(lin(i),min(i),nin(i),gama,v(1),d)                                 
!$OMP PARALLEL PRIVATE ( j, d ) SHARED ( xyz,e,v,gama,lin,min,nin )
!$OMP DO
         do j=1,nuc
            d(1:3) = e(1:3) - xyz(1:3,j)
            call opaa0(lin(i),min(i),nin(i),gama,v(j+1),d)                                 
         enddo
!$OMP END DO
!$OMP END PARALLEL      
         do j=1,nuc+1                                                          
            va(j)=va(j)+dd(i)*v(j)*efact                                                   
         enddo
      enddo
      return                                                                    
  120 continue                                                                  
      if(iff1.ne.1.and.iff2.ne.1) goto 121                                      
c           s - d                                                               
      do i=1,10                                                              
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         d = e 
         call opad1(lin(i),min(i),nin(i),gama,v(1),d)                                 
!$OMP PARALLEL PRIVATE ( j, d ) SHARED ( xyz,e,v,gama,lin,min,nin )
!$OMP DO
         do j=1,nuc
            d(1:3) = e(1:3) - xyz(1:3,j)
            call opaa0(lin(i),min(i),nin(i),gama,v(j+1),d)                                 
         enddo
!$OMP END DO
!$OMP END PARALLEL      
         do j=1,nuc+1                                                          
            va(j)=va(j)+dd(i)*v(j)*efact                                                   
         enddo
      enddo
      return                                                                    
  121 continue                                                                  
      if(iff1.gt.4.and.iff2.gt.4) goto 122                                      
c           p - d                                                               
      do i=1,20                                                              
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         d = e 
         call opad1(lin(i),min(i),nin(i),gama,v(1),d)                                 
!$OMP PARALLEL PRIVATE ( j, d ) SHARED ( xyz,e,v,gama,lin,min,nin )
!$OMP DO
         do j=1,nuc
            d(1:3) = e(1:3) - xyz(1:3,j)
            call opaa0(lin(i),min(i),nin(i),gama,v(j+1),d)                                 
         enddo
!$OMP END DO     
!$OMP END PARALLEL      
         do j=1,nuc+1                                                          
            va(j)=va(j)+dd(i)*v(j)*efact                                                   
         enddo
      enddo
      return                                                                    
  122 continue                                                                  
c           d - d                                                               
      do i=1,35                                                              
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         d = e 
         call opad1(lin(i),min(i),nin(i),gama,v(1),d)                                 
!$OMP PARALLEL PRIVATE ( j, d ) SHARED ( xyz,e,v,gama,lin,min,nin )
!$OMP DO
         do j=1,nuc
            d(1:3) = e(1:3) - xyz(1:3,j)
            call opaa0(lin(i),min(i),nin(i),gama,v(j+1),d)                                 
         enddo
!$OMP END DO     
!$OMP END PARALLEL      
         do j=1,nuc+1                                                          
            va(j)=va(j)+dd(i)*v(j)*efact                                                   
         enddo
      enddo
      return                                                                    
  110 continue                                                                  
      stop 'no f-fkt.'

      end                                                                       

c=======================================================================
c cartesian gaussian functions (6d,10f...)
c iff :
c s,px, py pz, dx**2 dy**2 dz**2 dxy dxz dyz
c 1 2   3   4   5     6     7     8   9  10
c fxxx, fyyy, fzzz, fxxy, fxxz, fyyx, fyyz, fxzz, fyzz, fxyz
c   11   12    13    14    15    16    17    18   19    20
c a, are aufpunkte
c nt : # of returns
c va : integral
c=======================================================================

      subroutine propa_st(nuc,xyz,a,b,etaij4,etakl4,iff1,iff2,va)   
      implicit real*8(a-h,o-z)                                                  
c aufpunkte,ref point,intarray
      real*8 a(3),b(3),va(2),xyz(3,nuc)                                     
c local
      common /abfunc/ ra(3),rb(3),ga,gb,ia,ib
      dimension d(3),dd(84),v(nuc+2)
      dimension e(3),aa(20),bb(20)
      dimension lin(84),min(84),nin(84)                                         
      data lin/0,1,0,0,2,0,0,1,1,0,3,0,0,2,2,1,0,1,0,1,4,0,0,3,3,1,0,1,0        
     1 ,2,2,0,2,1,1,5,0,0,3,3,2,2,0,0,4,4,1,0,0,1,1,3,1,2,2,1,6,0,0,3,          
     2 3,0,5,5,1,0,0,1,4,4,2,0,2,0,3,3,1,2,2,1,4,1,1,2/,min/0,0,1,0,0,          
     3 2,0,1,0,1,0,3,0,1,0,2,2,0,1,1,0,4,0,1,0,3,3,0,1,2,0,2,1,2,1,0,5,         
     4 0,2,0,3,0,3,2,1,0,4,4,1,0,1,1,3,2,1,2,0,6,0,3,0,3,1,0,0,1,5,5,2,0        
     5,0,2,4,4,2,1,3,1,3,2,1,4,1,2/,nin/0,0,0,1,0,0,2,0,1,1,0,0,3,0,1,0,        
     6 1,2,2,1,0,0,4,0,1,0,1,3,3,0,2,2,1,1,2,0,0,5,0,2,0,3,2,3,0,1,0,1,         
     7 4,4,3,1,1,1,2,2,0,0,6,0,3,3,0,1,5,5,1,0,0,2,4,4,0,2,1,2,2,3,1,3,         
     8 1,1,4,2/                                                                 

c --- a,b are centres of gaussians                                              
      ra = a
      rb = b
c --- ga,gb are their exponents                                                 
      ga=etaij4                                                              
      gb=etakl4                                                              
      aa = 0
      bb = 0
      aa(iff1)=1.0d0
      bb(iff2)=1.0d0
      ia=iff1                                                                   
      ib=iff2                                                                   
c --- ia,ib are the canonical indices for monom                                 
c --- apply product theorem                                                     
      cij=0.0d0
      ckl=0.0d0
      call divpt(a,etaij4,b,etakl4,cij,ckl,e,gama,efact)                  

c --- calculate cartesian prefactor for first gaussian                          
      call rhftce(aa,a,e,iff1)                                                  
c --- calculate cartesian prefactor for second gaussian                         
      call rhftce(bb,b,e,iff2)                                                  
c --- form their product                                                        
      call prod(aa,bb,dd,iff1,iff2)                                             
      va = 0
c ----- e is center of product gaussian with exponent gama                      
c          aname represents an external function                                
      if(iff1.gt.10.or.iff2.gt.10) goto 110                                     
      if(iff1.gt. 4.or.iff2.gt. 4) goto 120                                       
      if(iff1.gt. 1.or.iff2.gt. 1) goto 130                                       
c           s - s                                                               
      d = e 
      call opad1(lin(1),min(1),nin(1),gama,v(1),d)                                 
      call opap4(lin(1),min(1),nin(1),gama,v(2),d)                                 
      do j=1,2                                                          
         va(j)=dd(1)*v(j)*efact                                                   
      enddo
      return                                                                    
  130 continue                                                                  
      if(iff1.ne.1.and.iff2.ne.1) goto 131                                      
c           s - p                                                               
      do i=1,4                                                               
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         d = e 
         call opad1(lin(i),min(i),nin(i),gama,v(1),d)                                 
         call opap4(lin(i),min(i),nin(i),gama,v(2),d)                                 
         do j=1,2                                                          
!           va(j)=va(j)+dd(i)*v(j)*efact                                                   
            va(j)=va(j)+dd(i)*v(j)
         enddo
      enddo
      va(1:2)=va(1:2)*efact
      return                                                                    
  131 continue                                                                  
c           p - p                                                               
      do i=1,10                                                              
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         d = e 
         call opad1(lin(i),min(i),nin(i),gama,v(1),d)                                 
         call opap4(lin(i),min(i),nin(i),gama,v(2),d)                                 
         do j=1,2                                                          
!           va(j)=va(j)+dd(i)*v(j)*efact                                                   
            va(j)=va(j)+dd(i)*v(j)
         enddo
      enddo
      va(1:2)=va(1:2)*efact
      return                                                                    
  120 continue                                                                  
      if(iff1.ne.1.and.iff2.ne.1) goto 121                                      
c           s - d                                                               
      do i=1,10                                                              
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         d = e 
         call opad1(lin(i),min(i),nin(i),gama,v(1),d)                                 
         call opap4(lin(i),min(i),nin(i),gama,v(2),d)                                 
         do j=1,2                                                          
!           va(j)=va(j)+dd(i)*v(j)*efact                                                   
            va(j)=va(j)+dd(i)*v(j)
         enddo
      enddo
      va(1:2)=va(1:2)*efact
      return                                                                    
  121 continue                                                                  
      if(iff1.gt.4.and.iff2.gt.4) goto 122                                      
c           p - d                                                               
      do i=1,20                                                              
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         d = e 
         call opad1(lin(i),min(i),nin(i),gama,v(1),d)                                 
         call opap4(lin(i),min(i),nin(i),gama,v(2),d)                                 
         do j=1,2                                                          
!           va(j)=va(j)+dd(i)*v(j)*efact                                                   
            va(j)=va(j)+dd(i)*v(j)
         enddo
      enddo
      va(1:2)=va(1:2)*efact
      return                                                                    
  122 continue                                                                  
c           d - d                                                               
      do i=1,35                                                              
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         d = e 
         call opad1(lin(i),min(i),nin(i),gama,v(1),d)                                 
         call opap4(lin(i),min(i),nin(i),gama,v(2),d)                                 
         do j=1,2                                                          
!           va(j)=va(j)+dd(i)*v(j)*efact                                                   
            va(j)=va(j)+dd(i)*v(j)                                                         
         enddo
      enddo
      va(1:2)=va(1:2)*efact
      return                                                                    
  110 continue                                                                  
      stop 'no f-fkt.'

      end                                                                       

c=======================================================================
c cartesian gaussian functions (6d,10f...)
c iff :
c s,px, py pz, dx**2 dy**2 dz**2 dxy dxz dyz
c 1 2   3   4   5     6     7     8   9  10
c fxxx, fyyy, fzzz, fxxy, fxxz, fyyx, fyyz, fxzz, fyzz, fxyz
c   11   12    13    14    15    16    17    18   19    20
c a, are aufpunkte
c nt : # of returns
c va : integral
c=======================================================================

      subroutine propa_s(nuc,xyz,a,b,etaij4,etakl4,iff1,iff2,va)   
      implicit real*8(a-h,o-z)                                                  
c aufpunkte,ref point,intarray
      real*8 a(3),b(3),va(1),xyz(3,nuc)                                     
c local
      common /abfunc/ ra(3),rb(3),ga,gb,ia,ib
      dimension d(3),dd(84),v(nuc+2)
      dimension e(3),aa(20),bb(20)
      dimension lin(84),min(84),nin(84)                                         
      data lin/0,1,0,0,2,0,0,1,1,0,3,0,0,2,2,1,0,1,0,1,4,0,0,3,3,1,0,1,0        
     1 ,2,2,0,2,1,1,5,0,0,3,3,2,2,0,0,4,4,1,0,0,1,1,3,1,2,2,1,6,0,0,3,          
     2 3,0,5,5,1,0,0,1,4,4,2,0,2,0,3,3,1,2,2,1,4,1,1,2/,min/0,0,1,0,0,          
     3 2,0,1,0,1,0,3,0,1,0,2,2,0,1,1,0,4,0,1,0,3,3,0,1,2,0,2,1,2,1,0,5,         
     4 0,2,0,3,0,3,2,1,0,4,4,1,0,1,1,3,2,1,2,0,6,0,3,0,3,1,0,0,1,5,5,2,0        
     5,0,2,4,4,2,1,3,1,3,2,1,4,1,2/,nin/0,0,0,1,0,0,2,0,1,1,0,0,3,0,1,0,        
     6 1,2,2,1,0,0,4,0,1,0,1,3,3,0,2,2,1,1,2,0,0,5,0,2,0,3,2,3,0,1,0,1,         
     7 4,4,3,1,1,1,2,2,0,0,6,0,3,3,0,1,5,5,1,0,0,2,4,4,0,2,1,2,2,3,1,3,         
     8 1,1,4,2/                                                                 

c --- a,b are centres of gaussians                                              
      ra = a
      rb = b
c --- ga,gb are their exponents                                                 
      ga=etaij4                                                              
      gb=etakl4                                                              
      aa = 0
      bb = 0
      aa(iff1)=1.0d0
      bb(iff2)=1.0d0
      ia=iff1                                                                   
      ib=iff2                                                                   
c --- ia,ib are the canonical indices for monom                                 
c --- apply product theorem                                                     
      cij=0.0d0
      ckl=0.0d0
      call divpt(a,etaij4,b,etakl4,cij,ckl,e,gama,efact)                  

c --- calculate cartesian prefactor for first gaussian                          
      call rhftce(aa,a,e,iff1)                                                  
c --- calculate cartesian prefactor for second gaussian                         
      call rhftce(bb,b,e,iff2)                                                  
c --- form their product                                                        
      call prod(aa,bb,dd,iff1,iff2)                                             
      va = 0
c ----- e is center of product gaussian with exponent gama                      
c          aname represents an external function                                
      if(iff1.gt.10.or.iff2.gt.10) goto 110                                     
      if(iff1.gt. 4.or.iff2.gt. 4) goto 120                                       
      if(iff1.gt. 1.or.iff2.gt. 1) goto 130                                       
c           s - s                                                               
      d = e 
      call opad1(lin(1),min(1),nin(1),gama,v(1),d)                                 
      do j=1,1                                                          
         va(j)=dd(1)*v(j)*efact                                                   
      enddo
      return                                                                    
  130 continue                                                                  
      if(iff1.ne.1.and.iff2.ne.1) goto 131                                      
c           s - p                                                               
      do i=1,4                                                               
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         d = e 
         call opad1(lin(i),min(i),nin(i),gama,v(1),d)                                 
         do j=1,1                                                          
            va(j)=va(j)+dd(i)*v(j)
         enddo
      enddo
      va(1)=va(1)*efact
      return                                                                    
  131 continue                                                                  
c           p - p                                                               
      do i=1,10                                                              
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         d = e 
         call opad1(lin(i),min(i),nin(i),gama,v(1),d)                                 
         do j=1,1                                                          
            va(j)=va(j)+dd(i)*v(j)
         enddo
      enddo
      va(1)=va(1)*efact
      return                                                                    
  120 continue                                                                  
      if(iff1.ne.1.and.iff2.ne.1) goto 121                                      
c           s - d                                                               
      do i=1,10                                                              
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         d = e 
         call opad1(lin(i),min(i),nin(i),gama,v(1),d)                                 
         do j=1,1                                                          
            va(j)=va(j)+dd(i)*v(j)
         enddo
      enddo
      va(1)=va(1)*efact
      return                                                                    
  121 continue                                                                  
      if(iff1.gt.4.and.iff2.gt.4) goto 122                                      
c           p - d                                                               
      do i=1,20                                                              
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         d = e 
         call opad1(lin(i),min(i),nin(i),gama,v(1),d)                                 
         do j=1,1                                                          
            va(j)=va(j)+dd(i)*v(j)
         enddo
      enddo
      va(1)=va(1)*efact
      return                                                                    
  122 continue                                                                  
c           d - d                                                               
      do i=1,35                                                              
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         d = e 
         call opad1(lin(i),min(i),nin(i),gama,v(1),d)                                 
         do j=1,1                                                          
            va(j)=va(j)+dd(i)*v(j)
         enddo
      enddo
      va(1)=va(1)*efact
      return                                                                    
  110 continue                                                                  
      stop 'no f-fkt.'

      end                                                                       

c=======================================================================
c cartesian gaussian functions (6d,10f...)
c iff :
c s,px, py pz, dx**2 dy**2 dz**2 dxy dxz dyz
c 1 2   3   4   5     6     7     8   9  10
c fxxx, fyyy, fzzz, fxxy, fxxz, fyyx, fyyz, fxzz, fyzz, fxyz
c   11   12    13    14    15    16    17    18   19    20
c a, are aufpunkte
c nt : # of returns
c va : integral
c=======================================================================

      subroutine propa_dip(c,a,b,etaij4,etakl4,iff1,iff2,va)   
      implicit real*8(a-h,o-z)                                                  
c aufpunkte,ref point,intarray
      real*8 c(3),a(3),b(3),va(3)
c local
      common /abfunc/ ra(3),rb(3),ga,gb,ia,ib
      dimension d(3),dd(84),v(3)
      dimension aa(20),bb(20),e(3)
      dimension lin(84),min(84),nin(84)                                         
      data lin/0,1,0,0,2,0,0,1,1,0,3,0,0,2,2,1,0,1,0,1,4,0,0,3,3,1,0,1,0        
     1 ,2,2,0,2,1,1,5,0,0,3,3,2,2,0,0,4,4,1,0,0,1,1,3,1,2,2,1,6,0,0,3,          
     2 3,0,5,5,1,0,0,1,4,4,2,0,2,0,3,3,1,2,2,1,4,1,1,2/,min/0,0,1,0,0,          
     3 2,0,1,0,1,0,3,0,1,0,2,2,0,1,1,0,4,0,1,0,3,3,0,1,2,0,2,1,2,1,0,5,         
     4 0,2,0,3,0,3,2,1,0,4,4,1,0,1,1,3,2,1,2,0,6,0,3,0,3,1,0,0,1,5,5,2,0        
     5,0,2,4,4,2,1,3,1,3,2,1,4,1,2/,nin/0,0,0,1,0,0,2,0,1,1,0,0,3,0,1,0,        
     6 1,2,2,1,0,0,4,0,1,0,1,3,3,0,2,2,1,1,2,0,0,5,0,2,0,3,2,3,0,1,0,1,         
     7 4,4,3,1,1,1,2,2,0,0,6,0,3,3,0,1,5,5,1,0,0,2,4,4,0,2,1,2,2,3,1,3,         
     8 1,1,4,2/                                                                 

c --- a,b are centres of gaussians                                              
      ra = a
      rb = b
c --- ga,gb are their exponents                                                 
      ga=etaij4                                                              
      gb=etakl4                                                              
      aa = 0
      bb = 0
      aa(iff1)=1.0d0
      bb(iff2)=1.0d0
      ia=iff1                                                                   
      ib=iff2                                                                   
c --- ia,ib are the canonical indices for monom                                 
c --- apply product theorem                                                     
      cij=0.0d0
      ckl=0.0d0
      call divpt(a,etaij4,b,etakl4,cij,ckl,e,gama,efact)                  

c --- calculate cartesian prefactor for first gaussian                          
      call rhftce(aa,a,e,iff1)                                                  
c --- calculate cartesian prefactor for second gaussian                         
      call rhftce(bb,b,e,iff2)                                                  
c --- form their product                                                        
      call prod(aa,bb,dd,iff1,iff2)                                             
      va = 0
c ----- e is center of product gaussian with exponent gama                      
c ---- reference point is c
      d  = e - c
c          aname represents an external function                                
      if(iff1.gt.10.or.iff2.gt.10) goto 110                                     
      if(iff1.gt. 4.or.iff2.gt. 4) goto 120                                       
      if(iff1.gt. 1.or.iff2.gt. 1) goto 130                                       
c           s - s                                                               
      call opab1(lin(1),min(1),nin(1),gama,v,d)                                 
      do j=1,3                                                          
         va(j)=-dd(1)*v(j)*efact                                                   
      enddo
      return                                                                    
  130 continue                                                                  
      if(iff1.ne.1.and.iff2.ne.1) goto 131                                      
c           s - p                                                               
      do i=1,4                                                               
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         call opab1(lin(i),min(i),nin(i),gama,v,d)                                 
         do j=1,3                                                          
            va(j)=va(j)-dd(i)*v(j)
         enddo
      enddo
      va(1:3)=va(1:3)*efact
      return                                                                    
  131 continue                                                                  
c           p - p                                                               
      do i=1,10                                                              
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         call opab1(lin(i),min(i),nin(i),gama,v,d)                                 
         do j=1,3                                                          
            va(j)=va(j)-dd(i)*v(j)
         enddo
      enddo
      va(1:3)=va(1:3)*efact
      return                                                                    
  120 continue                                                                  
      if(iff1.ne.1.and.iff2.ne.1) goto 121                                      
c           s - d                                                               
      do i=1,10                                                              
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         call opab1(lin(i),min(i),nin(i),gama,v,d)                                 
         do j=1,3                                                          
            va(j)=va(j)-dd(i)*v(j)
         enddo
      enddo
      va(1:3)=va(1:3)*efact
      return                                                                    
  121 continue                                                                  
      if(iff1.gt.4.and.iff2.gt.4) goto 122                                      
c           p - d                                                               
      do i=1,20                                                              
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         call opab1(lin(i),min(i),nin(i),gama,v,d)                                 
         do j=1,3                                                          
            va(j)=va(j)-dd(i)*v(j)
         enddo
      enddo
      va(1:3)=va(1:3)*efact
      return                                                                    
  122 continue                                                                  
c           d - d                                                               
      do i=1,35                                                              
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         call opab1(lin(i),min(i),nin(i),gama,v,d)                                 
         do j=1,3                                                          
            va(j)=va(j)-dd(i)*v(j)
         enddo
      enddo
      va(1:3)=va(1:3)*efact
      return                                                                    
  110 continue                                                                  
      stop 'no f-fkt.'

      end                                                                       

c=======================================================================
c x^2,y^2,z^2 only
c=======================================================================

      subroutine propa_sec(c,a,b,etaij4,etakl4,iff1,iff2,va)   
      implicit real*8(a-h,o-z)                                                  
c aufpunkte,ref point,intarray
      real*8 c(3),a(3),b(3),va(3)
c local
      common /abfunc/ ra(3),rb(3),ga,gb,ia,ib
      dimension d(3),dd(84),v(3)
      dimension aa(20),bb(20),e(3)
      dimension lin(84),min(84),nin(84)                                         
      data lin/0,1,0,0,2,0,0,1,1,0,3,0,0,2,2,1,0,1,0,1,4,0,0,3,3,1,0,1,0        
     1 ,2,2,0,2,1,1,5,0,0,3,3,2,2,0,0,4,4,1,0,0,1,1,3,1,2,2,1,6,0,0,3,          
     2 3,0,5,5,1,0,0,1,4,4,2,0,2,0,3,3,1,2,2,1,4,1,1,2/,min/0,0,1,0,0,          
     3 2,0,1,0,1,0,3,0,1,0,2,2,0,1,1,0,4,0,1,0,3,3,0,1,2,0,2,1,2,1,0,5,         
     4 0,2,0,3,0,3,2,1,0,4,4,1,0,1,1,3,2,1,2,0,6,0,3,0,3,1,0,0,1,5,5,2,0        
     5,0,2,4,4,2,1,3,1,3,2,1,4,1,2/,nin/0,0,0,1,0,0,2,0,1,1,0,0,3,0,1,0,        
     6 1,2,2,1,0,0,4,0,1,0,1,3,3,0,2,2,1,1,2,0,0,5,0,2,0,3,2,3,0,1,0,1,         
     7 4,4,3,1,1,1,2,2,0,0,6,0,3,3,0,1,5,5,1,0,0,2,4,4,0,2,1,2,2,3,1,3,         
     8 1,1,4,2/                                                                 

c --- a,b are centres of gaussians                                              
      ra = a
      rb = b
c --- ga,gb are their exponents                                                 
      ga=etaij4                                                              
      gb=etakl4                                                              
      aa = 0
      bb = 0
      aa(iff1)=1.0d0
      bb(iff2)=1.0d0
      ia=iff1                                                                   
      ib=iff2                                                                   
c --- ia,ib are the canonical indices for monom                                 
c --- apply product theorem                                                     
      cij=0.0d0
      ckl=0.0d0
      call divpt(a,etaij4,b,etakl4,cij,ckl,e,gama,efact)                  

c --- calculate cartesian prefactor for first gaussian                          
      call rhftce(aa,a,e,iff1)                                                  
c --- calculate cartesian prefactor for second gaussian                         
      call rhftce(bb,b,e,iff2)                                                  
c --- form their product                                                        
      call prod(aa,bb,dd,iff1,iff2)                                             
      va = 0
c ----- e is center of product gaussian with exponent gama                      
c ---- reference point is c
      d  = e - c
c          aname represents an external function                                
      if(iff1.gt.10.or.iff2.gt.10) goto 110                                     
      if(iff1.gt. 4.or.iff2.gt. 4) goto 120                                       
      if(iff1.gt. 1.or.iff2.gt. 1) goto 130                                       
c           s - s                                                               
      call opab4(lin(1),min(1),nin(1),gama,v,d)                                 
      do j=1,3                                                          
         va(j)=dd(1)*v(j)*efact                                                   
      enddo
      return                                                                    
  130 continue                                                                  
      if(iff1.ne.1.and.iff2.ne.1) goto 131                                      
c           s - p                                                               
      do i=1,4                                                               
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         call opab4(lin(i),min(i),nin(i),gama,v,d)                                 
         do j=1,3                                                          
            va(j)=va(j)+dd(i)*v(j)
         enddo
      enddo
      va(1:3)=va(1:3)*efact
      return                                                                    
  131 continue                                                                  
c           p - p                                                               
      do i=1,10                                                              
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         call opab4(lin(i),min(i),nin(i),gama,v,d)                                 
         do j=1,3                                                          
            va(j)=va(j)+dd(i)*v(j)
         enddo
      enddo
      va(1:3)=va(1:3)*efact
      return                                                                    
  120 continue                                                                  
      if(iff1.ne.1.and.iff2.ne.1) goto 121                                      
c           s - d                                                               
      do i=1,10                                                              
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         call opab4(lin(i),min(i),nin(i),gama,v,d)                                 
         do j=1,3                                                          
            va(j)=va(j)+dd(i)*v(j)
         enddo
      enddo
      va(1:3)=va(1:3)*efact
      return                                                                    
  121 continue                                                                  
      if(iff1.gt.4.and.iff2.gt.4) goto 122                                      
c           p - d                                                               
      do i=1,20                                                              
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         call opab4(lin(i),min(i),nin(i),gama,v,d)                                 
         do j=1,3                                                          
            va(j)=va(j)+dd(i)*v(j)
         enddo
      enddo
      va(1:3)=va(1:3)*efact
      return                                                                    
  122 continue                                                                  
c           d - d                                                               
      do i=1,35                                                              
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         call opab4(lin(i),min(i),nin(i),gama,v,d)                                 
         do j=1,3                                                          
            va(j)=va(j)+dd(i)*v(j)
         enddo
      enddo
      va(1:3)=va(1:3)*efact
      return                                                                    
  110 continue                                                                  
      stop 'no f-fkt.'

      end                                                                       

      subroutine propa_dq(c,a,b,etaij4,etakl4,iff1,iff2,va)   
      implicit real*8(a-h,o-z)                                                  
c aufpunkte,ref point,intarray
      real*8 c(3),a(3),b(3),va(9)
c local
      common /abfunc/ ra(3),rb(3),ga,gb,ia,ib
      dimension d(3),dd(84),v(9)
      dimension aa(20),bb(20),e(3)
      dimension lin(84),min(84),nin(84)                                         
      data lin/0,1,0,0,2,0,0,1,1,0,3,0,0,2,2,1,0,1,0,1,4,0,0,3,3,1,0,1,0        
     1 ,2,2,0,2,1,1,5,0,0,3,3,2,2,0,0,4,4,1,0,0,1,1,3,1,2,2,1,6,0,0,3,          
     2 3,0,5,5,1,0,0,1,4,4,2,0,2,0,3,3,1,2,2,1,4,1,1,2/,min/0,0,1,0,0,          
     3 2,0,1,0,1,0,3,0,1,0,2,2,0,1,1,0,4,0,1,0,3,3,0,1,2,0,2,1,2,1,0,5,         
     4 0,2,0,3,0,3,2,1,0,4,4,1,0,1,1,3,2,1,2,0,6,0,3,0,3,1,0,0,1,5,5,2,0        
     5,0,2,4,4,2,1,3,1,3,2,1,4,1,2/,nin/0,0,0,1,0,0,2,0,1,1,0,0,3,0,1,0,        
     6 1,2,2,1,0,0,4,0,1,0,1,3,3,0,2,2,1,1,2,0,0,5,0,2,0,3,2,3,0,1,0,1,         
     7 4,4,3,1,1,1,2,2,0,0,6,0,3,3,0,1,5,5,1,0,0,2,4,4,0,2,1,2,2,3,1,3,         
     8 1,1,4,2/                                                                 

c --- a,b are centres of gaussians                                              
      ra = a
      rb = b
c --- ga,gb are their exponents                                                 
      ga=etaij4                                                              
      gb=etakl4                                                              
      aa = 0
      bb = 0
      aa(iff1)=1.0d0
      bb(iff2)=1.0d0
      ia=iff1                                                                   
      ib=iff2                                                                   
c --- ia,ib are the canonical indices for monom                                 
c --- apply product theorem                                                     
      cij=0.0d0
      ckl=0.0d0
      call divpt(a,etaij4,b,etakl4,cij,ckl,e,gama,efact)                  

c --- calculate cartesian prefactor for first gaussian                          
      call rhftce(aa,a,e,iff1)                                                  
c --- calculate cartesian prefactor for second gaussian                         
      call rhftce(bb,b,e,iff2)                                                  
c --- form their product                                                        
      call prod(aa,bb,dd,iff1,iff2)                                             
      va = 0
c ----- e is center of product gaussian with exponent gama                      
c ---- reference point is c
      d  = e - c
c          aname represents an external function                                
      if(iff1.gt.10.or.iff2.gt.10) goto 110                                     
      if(iff1.gt. 4.or.iff2.gt. 4) goto 120                                       
      if(iff1.gt. 1.or.iff2.gt. 1) goto 130                                       
c           s - s                                                               
      call opab1 (lin(1),min(1),nin(1),gama,v,d)                                 
      call opab46(lin(1),min(1),nin(1),gama,v,d)                                 
      do j=1,9                                                          
         va(j)=dd(1)*v(j)*efact                                                   
      enddo
      return                                                                    
  130 continue                                                                  
      if(iff1.ne.1.and.iff2.ne.1) goto 131                                      
c           s - p                                                               
      do i=1,4                                                               
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         call opab1 (lin(i),min(i),nin(i),gama,v,d)                                 
         call opab46(lin(i),min(i),nin(i),gama,v,d)                                 
         do j=1,9                                                          
            va(j)=va(j)+dd(i)*v(j)
         enddo
      enddo
      va(1:9)=va(1:9)*efact
      return                                                                    
  131 continue                                                                  
c           p - p                                                               
      do i=1,10                                                              
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         call opab1 (lin(i),min(i),nin(i),gama,v,d)                                 
         call opab46(lin(i),min(i),nin(i),gama,v,d)                                 
         do j=1,9                                                          
            va(j)=va(j)+dd(i)*v(j)
         enddo
      enddo
      va(1:9)=va(1:9)*efact
      return                                                                    
  120 continue                                                                  
      if(iff1.ne.1.and.iff2.ne.1) goto 121                                      
c           s - d                                                               
      do i=1,10                                                              
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         call opab1 (lin(i),min(i),nin(i),gama,v,d)                                 
         call opab46(lin(i),min(i),nin(i),gama,v,d)                                 
         do j=1,9                                                          
            va(j)=va(j)+dd(i)*v(j)
         enddo
      enddo
      va(1:9)=va(1:9)*efact
      return                                                                    
  121 continue                                                                  
      if(iff1.gt.4.and.iff2.gt.4) goto 122                                      
c           p - d                                                               
      do i=1,20                                                              
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         call opab1 (lin(i),min(i),nin(i),gama,v,d)                                 
         call opab46(lin(i),min(i),nin(i),gama,v,d)                                 
         do j=1,9                                                          
            va(j)=va(j)+dd(i)*v(j)
         enddo
      enddo
      va(1:9)=va(1:9)*efact
      return                                                                    
  122 continue                                                                  
c           d - d                                                               
      do i=1,35                                                              
         if(dabs(dd(i))-1.d-8.le.0) cycle                                               
         call opab1 (lin(i),min(i),nin(i),gama,v,d)                                 
         call opab46(lin(i),min(i),nin(i),gama,v,d)                                 
         do j=1,9                                                          
            va(j)=va(j)+dd(i)*v(j)
         enddo
      enddo
      va(1:9)=va(1:9)*efact
      return                                                                    
  110 continue                                                                  
      stop 'no f-fkt.'

      end                                                                       

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine divpt(a,alpha,b,beta,ca,cb,e,fgama,ffact)                      
      implicit real*8(a-h,o-z)                                                  
c         this computes the center, exponent, and multiplying factor of         
c         a single gaussian which can replace the product of two gaussia        
c         centers a and b, and exponents alpha and beta.                        
      dimension a(3), b(3), e(3)                                                
      data  toluol /20.7232660d0/                                                     
      abexp=alpha+beta                                                          
      do 1 i=1,3                                                                
    1 e(i)=(alpha*a(i)+beta*b(i))/abexp                                         
      absqd=0.0d0                                                               
      do 2 i=1,3                                                                
    2 absqd=absqd+(b(i)-a(i))*(b(i)-a(i))                                       
      abprod=absqd*alpha*beta/abexp                                             
      fgama=abexp                                                               
      ffact=0.d0                                                                
      tol=toluol+ca+cb                                                          
      if(abprod.gt.tol) return                                                  
      ffact=dexp(-abprod)                                                       
      return                                                                    
      end                                                                       

      subroutine rhftce(cfs,a,e,iff)                                            
      implicit real*8(a-h,o-z)                                                  
      dimension cfs(*) ,a(*),e(*)                                               
      data c2/2.0d0/,c3/3.0d0/                                                  
c ---- e = center of product function, a = center of single gaussian            
      aex = e(1)-a(1)                                                           
      aey = e(2)-a(2)                                                           
      aez = e(3)-a(3)                                                           
      goto (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),iff             
    1 return                                                                    
    2 cfs(1)=aex*cfs(2)                                                         
      return                                                                    
    3 cfs(1)=aey*cfs(3)                                                         
      return                                                                    
    4 cfs(1)=aez*cfs(4)                                                         
      return                                                                    
    5 cfs(1)=aex*aex*cfs(5)                                                     
      cfs(2)=c2*aex*cfs(5)                                                      
      return                                                                    
    6 cfs(1)=aey*aey*cfs(6)                                                     
      cfs(3)=c2*aey*cfs(6)                                                      
      return                                                                    
    7 cfs(1)=aez*aez*cfs(7)                                                     
      cfs(4)=c2*aez*cfs(7)                                                      
      return                                                                    
    8 cfs(1)=aex*aey*cfs(8)                                                     
      cfs(2)=aey*cfs(8)                                                         
      cfs(3)=aex*cfs(8)                                                         
      return                                                                    
    9 cfs(1)=aex*aez*cfs(9)                                                     
      cfs(2)=aez*cfs(9)                                                         
      cfs(4)=aex*cfs(9)                                                         
      return                                                                    
   10 cfs(1)=aey*aez*cfs(10)                                                    
      cfs(3)=aez*cfs(10)                                                        
      cfs(4)=aey*cfs(10)                                                        
      return                                                                    
   11 cfs(1)=aex*aex*aex*cfs(11)                                                
      cfs(2)=c3*aex*aex*cfs(11)                                                 
      cfs(5)=c3*aex*cfs(11)                                                     
      return                                                                    
   12 cfs(1)=aey*aey*aey*cfs(12)                                                
      cfs(3)=c3*aey*aey*cfs(12)                                                 
      cfs(6)=c3*aey*cfs(12)                                                     
      return                                                                    
   13 cfs(1)=aez*aez*aez*cfs(13)                                                
      cfs(4)=c3*aez*aez*cfs(13)                                                 
      cfs(7)=c3*aez*cfs(13)                                                     
      return                                                                    
   14 cfs(1)=aex*aex*aey*cfs(14)                                                
      cfs(2)=c2*aex*aey*cfs(14)                                                 
      cfs(3)=aex*aex*cfs(14)                                                    
      cfs(5)=aey*cfs(14)                                                        
      cfs(8)=c2*aex*cfs(14)                                                     
      return                                                                    
   15 cfs(1)=aex*aex*aez*cfs(15)                                                
      cfs(2)=c2*aex*aez*cfs(15)                                                 
      cfs(4)=aex*aex*cfs(15)                                                    
      cfs(5)=aez*cfs(15)                                                        
      cfs(9)=c2*aex*cfs(15)                                                     
      return                                                                    
   16 cfs(1)=aey*aey*aex*cfs(16)                                                
      cfs(2)=aey*aey*cfs(16)                                                    
      cfs(3)=c2*aey*aex*cfs(16)                                                 
      cfs(6)=aex*cfs(16)                                                        
      cfs(8)=c2*aey*cfs(16)                                                     
      return                                                                    
   17 cfs(1)=aey*aey*aez*cfs(17)                                                
      cfs(3)=c2*aey*aez*cfs(17)                                                 
      cfs(4)=aey*aey*cfs(17)                                                    
      cfs(6)=aez*cfs(17)                                                        
      cfs(10)=c2*aey*cfs(17)                                                    
      return                                                                    
   18 cfs(1)=aez*aez*aex*cfs(18)                                                
      cfs(2)=aez*aez*cfs(18)                                                    
      cfs(4)=c2*aez*aex*cfs(18)                                                 
      cfs(7)=aex*cfs(18)                                                        
      cfs(9)=c2*aez*cfs(18)                                                     
      return                                                                    
   19 cfs(1)=aez*aez*aey*cfs(19)                                                
      cfs(3)=aez*aez*cfs(19)                                                    
      cfs(4)=c2*aez*aey*cfs(19)                                                 
      cfs(7)=aey*cfs(19)                                                        
      cfs(10)=c2*aez*cfs(19)                                                    
      return                                                                    
   20 cfs(1)=aex*aey*aez*cfs(20)                                                
      cfs(2)=aez*aey*cfs(20)                                                    
      cfs(3)=aex*aez*cfs(20)                                                    
      cfs(4)=aex*aey*cfs(20)                                                    
      cfs(8)=aez*cfs(20)                                                        
      cfs(9)=aey*cfs(20)                                                        
      cfs(10)=aex*cfs(20)                                                       
      return                                                                    
      end                                                                       

      subroutine prod(c,d,s,iff1,iff2)                                          
      implicit real*8(a-h,o-z)                                                  
      dimension c(*),d(*),s(*)                                                  
      if(iff1.gt.10.or.iff2.gt.10) goto 30                                      
      if(iff1.gt.4.or.iff2.gt.4) goto 20                                        
      s( 1)=c( 1)*d( 1)                                                         
c              end of s - s                                                     
      if(iff1.eq.1.and.iff2.eq.1) return                                        
      s( 2)=c( 1)*d( 2)+c( 2)*d( 1)                                             
      s( 3)=c( 1)*d( 3)+c( 3)*d( 1)                                             
      s( 4)=c( 1)*d( 4)+c( 4)*d( 1)                                             
c              end of s - p                                                     
      if(iff1.eq.1.or.iff2.eq.1) return                                         
      s( 5)=c( 2)*d( 2)                                                         
      s( 6)=c( 3)*d( 3)                                                         
      s( 7)=c( 4)*d( 4)                                                         
      s( 8)=c( 2)*d( 3)+c( 3)*d( 2)                                             
      s( 9)=c( 2)*d( 4)+c( 4)*d( 2)                                             
      s(10)=c( 3)*d( 4)+c( 4)*d( 3)                                             
c              end of p - p                                                     
      return                                                                    
   20 continue                                                                  
      s( 1)=c( 1)*d( 1)                                                         
      s( 2)=c( 1)*d( 2)+c( 2)*d( 1)                                             
      s( 3)=c( 1)*d( 3)+c( 3)*d( 1)                                             
      s( 4)=c( 1)*d( 4)+c( 4)*d( 1)                                             
      s( 5)=c( 1)*d( 5)+c( 5)*d( 1)+c( 2)*d( 2)                                 
      s( 6)=c( 1)*d( 6)+c( 6)*d( 1)+c( 3)*d( 3)                                 
      s( 7)=c( 1)*d( 7)+c( 7)*d( 1)+c( 4)*d( 4)                                 
      s( 8)=c( 1)*d( 8)+c( 8)*d( 1)+c( 2)*d( 3)+c( 3)*d( 2)                     
      s( 9)=c( 1)*d( 9)+c( 9)*d( 1)+c( 2)*d( 4)+c( 4)*d( 2)                     
      s(10)=c( 1)*d(10)+c(10)*d( 1)+c( 3)*d( 4)+c( 4)*d( 3)                     
c              end of s - d                                                     
      if(iff1.eq.1.or.iff2.eq.1) return                                         
      s(11)=c( 2)*d( 5)+c( 5)*d( 2)                                             
      s(12)=c( 3)*d( 6)+c( 6)*d( 3)                                             
      s(13)=c( 4)*d( 7)+c( 7)*d( 4)                                             
      s(14)=c( 2)*d( 8)+c( 8)*d( 2)+c( 3)*d( 5)+c( 5)*d( 3)                     
      s(15)=c( 2)*d( 9)+c( 9)*d( 2)+c( 4)*d( 5)+c( 5)*d( 4)                     
      s(16)=c( 2)*d( 6)+c( 6)*d( 2)+c( 3)*d( 8)+c( 8)*d( 3)                     
      s(17)=c( 3)*d(10)+c(10)*d( 3)+c( 4)*d( 6)+c( 6)*d( 4)                     
      s(18)=c( 2)*d( 7)+c( 7)*d( 2)+c( 4)*d( 9)+c( 9)*d( 4)                     
      s(19)=c( 3)*d( 7)+c( 7)*d( 3)+c( 4)*d(10)+c(10)*d( 4)                     
      s(20)=c( 2)*d(10)+c(10)*d( 2)+c( 3)*d( 9)                                 
     1     +c( 9)*d( 3)+c( 4)*d( 8)+c( 8)*d( 4)                                 
c              end of p - d                                                     
      if(iff1.lt.5.or.iff2.lt.5) return                                         
      s(21)=c( 5)*d( 5)                                                         
      s(22)=c( 6)*d( 6)                                                         
      s(23)=c( 7)*d( 7)                                                         
      s(24)=c( 5)*d( 8)+c( 8)*d( 5)                                             
      s(25)=c( 5)*d( 9)+c( 9)*d( 5)                                             
      s(26)=c( 6)*d( 8)+c( 8)*d( 6)                                             
      s(27)=c( 6)*d(10)+c(10)*d( 6)                                             
      s(28)=c( 7)*d( 9)+c( 9)*d( 7)                                             
      s(29)=c( 7)*d(10)+c(10)*d( 7)                                             
      s(30)=c( 5)*d( 6)+c( 6)*d( 5)+c( 8)*d( 8)                                 
      s(31)=c( 5)*d( 7)+c( 7)*d( 5)+c( 9)*d( 9)                                 
      s(32)=c( 6)*d( 7)+c( 7)*d( 6)+c(10)*d(10)                                 
      s(33)=c( 5)*d(10)+c(10)*d( 5)+c( 8)*d( 9)+c( 9)*d( 8)                     
      s(34)=c( 6)*d( 9)+c( 9)*d( 6)+c( 8)*d(10)+c(10)*d( 8)                     
      s(35)=c( 7)*d( 8)+c( 8)*d( 7)+c( 9)*d(10)+c(10)*d( 9)                     
c              end of d - d                                                     
      return                                                                    
   30 continue                                                                  
      s( 1)=c( 1)*d( 1)                                                         
      s( 2)=c( 1)*d( 2)+c( 2)*d( 1)                                             
      s( 3)=c( 1)*d( 3)+c( 3)*d( 1)                                             
      s( 4)=c( 1)*d( 4)+c( 4)*d( 1)                                             
      s( 5)=c( 1)*d( 5)+c( 5)*d( 1)+c( 2)*d( 2)                                 
      s( 6)=c( 1)*d( 6)+c( 6)*d( 1)+c( 3)*d( 3)                                 
      s( 7)=c( 1)*d( 7)+c( 7)*d( 1)+c( 4)*d( 4)                                 
      s( 8)=c( 1)*d( 8)+c( 8)*d( 1)+c( 2)*d( 3)+c( 3)*d( 2)                     
      s( 9)=c( 1)*d( 9)+c( 9)*d( 1)+c( 2)*d( 4)+c( 4)*d( 2)                     
      s(10)=c( 1)*d(10)+c(10)*d( 1)+c( 3)*d( 4)+c( 4)*d( 3)                     
      s(11)=c( 2)*d( 5)+c( 5)*d( 2)+c(1)*d(11)+c(11)*d(1)                       
      s(12)=c( 3)*d( 6)+c( 6)*d( 3)+c(1)*d(12)+c(12)*d(1)                       
      s(13)=c( 4)*d( 7)+c( 7)*d( 4)+c(1)*d(13)+c(13)*d(1)                       
      s(14)=c( 2)*d( 8)+c( 8)*d( 2)+c( 3)*d( 5)+c( 5)*d( 3)+c(1)*d(14)+         
     1      c(14)*d(1)                                                          
      s(15)=c( 2)*d( 9)+c( 9)*d( 2)+c( 4)*d( 5)+c( 5)*d( 4)+c(1)*d(15)+         
     1      c(15)*d(1)                                                          
      s(16)=c( 2)*d( 6)+c( 6)*d( 2)+c( 3)*d( 8)+c( 8)*d( 3)+c(1)*d(16)+         
     1      c(16)*d(1)                                                          
      s(17)=c( 3)*d(10)+c(10)*d( 3)+c( 4)*d( 6)+c( 6)*d( 4)+c(1)*d(17)+         
     1      c(17)*d(1)                                                          
      s(18)=c( 2)*d( 7)+c( 7)*d( 2)+c( 4)*d( 9)+c( 9)*d( 4)+c(1)*d(18)+         
     1      c(18)*d(1)                                                          
      s(19)=c( 3)*d( 7)+c( 7)*d( 3)+c( 4)*d(10)+c(10)*d( 4)+c(1)*d(19)+         
     1      c(19)*d(1)                                                          
      s(20)=c( 2)*d(10)+c(10)*d( 2)+c( 3)*d( 9)+c(9)*d(3)+c(4)*d(8)+            
     1      c(8)*d(4)+c(1)*d(20)+c(20)*d(1)                                     
c              end of s - f                                                     
      if(iff1.eq.1.or.iff2.eq.1) return                                         
      s(21)=c( 5)*d( 5)+c(2)*d(11)+c(11)*d(2)                                   
      s(22)=c( 6)*d( 6)+c(3)*d(12)+c(12)*d(3)                                   
      s(23)=c( 7)*d( 7)+c(4)*d(13)+c(13)*d(4)                                   
      s(24)=c( 5)*d( 8)+c( 8)*d( 5)+c(3)*d(11)+c(11)*d(3)+c(2)*d(14)+           
     1      c(14)*d(2)                                                          
      s(25)=c( 5)*d( 9)+c( 9)*d( 5)+c(2)*d(15)+c(15)*d(2)+c(4)*d(11)+           
     1      c(11)*d(4)                                                          
      s(26)=c( 6)*d( 8)+c( 8)*d( 6)+c(2)*d(12)+c(12)*d(2)+c(3)*d(16)+           
     1      c(16)*d(3)                                                          
      s(27)=c( 6)*d(10)+c(10)*d( 6)+c(3)*d(17)+c(17)*d(3)+c(4)*d(12)+           
     1      c(12)*d(4)                                                          
      s(28)=c( 7)*d( 9)+c( 9)*d( 7)+c(2)*d(13)+c(13)*d(2)+c(4)*d(18)+           
     1      c(18)*d(4)                                                          
      s(29)=c( 7)*d(10)+c(10)*d( 7)+c(3)*d(13)+c(13)*d(3)+c(4)*d(19)+           
     1      c(19)*d(4)                                                          
      s(30)=c( 5)*d( 6)+c( 6)*d( 5)+c( 8)*d( 8)+c(2)*d(16)+c(16)*d(2)+          
     1      c(3)*d(14)+c(14)*d(3)                                               
      s(31)=c( 5)*d( 7)+c( 7)*d( 5)+c( 9)*d( 9)+c(2)*d(18)+c(18)*d(2)+          
     1      c(4)*d(15)+c(15)*d(4)                                               
      s(32)=c( 6)*d( 7)+c( 7)*d( 6)+c(10)*d(10)+c(3)*d(19)+c(19)*d(3)+          
     1      c(4)*d(17)+c(17)*d(4)                                               
      s(33)=c( 5)*d(10)+c(10)*d( 5)+c( 8)*d( 9)+c( 9)*d( 8)+c(3)*d(15)+         
     1      c(15)*d(3)+c(4)*d(14)+c(14)*d(4)+c(2)*d(20)+c(20)*d(2)              
      s(34)=c( 6)*d( 9)+c( 9)*d( 6)+c( 8)*d(10)+c(10)*d( 8)+c(2)*d(17)+         
     1      d(2)*c(17)+c(3)*d(20)+c(20)*d(3)+c(4)*d(16)+c(16)*d(4)              
      s(35)=c( 7)*d( 8)+c( 8)*d( 7)+c( 9)*d(10)+c(10)*d( 9)+c(2)*d(19)+         
     1      c(19)*d(2)+c(3)*d(18)+c(18)*d(3)+c(4)*d(20)+c(20)*d(4)              
c              end of p - f                                                     
      if(iff1.eq.2.or.iff2.eq.2) return                                         
      s(36)=c(5)*d(11)+c(11)*d(5)                                               
      s(37)=c(6)*d(12)+c(12)*d(6)                                               
      s(38)=c(7)*d(13)+c(13)*d(7)                                               
      s(39)=c(6)*d(11)+c(11)*d(6)+c(5)*d(16)+c(16)*d(5)+c(8)*d(14)+             
     1      c(14)*d(8)                                                          
      s(40)=c(7)*d(11)+c(11)*d(7)+c(5)*d(18)+c(18)*d(5)+c(9)*d(15)+             
     1      c(15)*d(9)                                                          
      s(41)=c(5)*d(12)+c(12)*d(5)+c(6)*d(14)+c(14)*d(6)+c(8)*d(16)+             
     1      c(16)*d(8)                                                          
      s(42)=c(5)*d(13)+c(13)*d(5)+c(7)*d(15)+c(15)*d(7)+c(9)*d(18)+             
     1      c(18)*d(9)                                                          
      s(43)=c(7)*d(12)+c(12)*d(7)+c(6)*d(19)+c(19)*d(6)+c(10)*d(17)+            
     1      c(17)*d(10)                                                         
      s(44)=c(6)*d(13)+c(13)*d(6)+c(7)*d(17)+c(17)*d(7)+c(10)*d(19)+            
     1      c(19)*d(10)                                                         
      s(45)=c(8)*d(11)+c(11)*d(8)+c(5)*d(14)+c(14)*d(5)                         
      s(46)=c(9)*d(11)+c(11)*d(9)+c(5)*d(15)+c(15)*d(5)                         
      s(47)=c(8)*d(12)+c(12)*d(8)+c(6)*d(16)+c(16)*d(6)                         
      s(48)=c(10)*d(12)+c(12)*d(10)+c(6)*d(17)+c(17)*d(6)                       
      s(49)=c(10)*d(13)+c(13)*d(10)+c(7)*d(19)+c(19)*d(7)                       
      s(50)=c(9)*d(13)+c(13)*d(9)+c(7)*d(18)+c(18)*d(7)                         
      s(51)=c(8)*d(13)+c(13)*d(8)+c(7)*d(20)+c(20)*d(7)+c(9)*d(19)+             
     1      c(19)*d(9)+c(10)*d(18)+c(18)*d(10)                                  
      s(52)=c(10)*d(11)+c(11)*d(10)+c(5)*d(20)+c(20)*d(5)+c(9)*d(14)+           
     1      c(14)*d(9)+c(8)*d(15)+c(15)*d(8)                                    
      s(53)=c(9)*d(12)+c(12)*d(9)+c(6)*d(20)+c(20)*d(6)+c(10)*d(16)+            
     1      c(16)*d(10)+c(8)*d(17)+c(17)*d(8)                                   
      s(54)=c(5)*d(17)+c(17)*d(5)+c(6)*d(15)+c(15)*d(6)+c(14)*d(10)+            
     1      d(14)*c(10)+c(9)*d(16)+c(16)*d(9)+c(8)*d(20)+c(20)*d(8)             
      s(55)=c(5)*d(19)+c(19)*d(5)+c(7)*d(14)+c(14)*d(7)+c(10)*d(15)+            
     1      c(15)*d(10)+c(9)*d(20)+c(20)*d(9)+c(8)*d(18)+c(18)*d(8)             
      s(56)=c(6)*d(18)+c(18)*d(6)+c(7)*d(16)+c(16)*d(7)+c(10)*d(20)+            
     1      c(20)*d(10)+c(9)*d(17)+c(17)*d(9)+c(8)*d(19)+c(19)*d(8)             
c              end of d - f                                                     
      if(iff1.eq.3.or.iff2.eq.3) return                                         
      s(57)=c(11)*d(11)                                                         
      s(58)=c(12)*d(12)                                                         
      s(59)=c(13)*d(13)                                                         
      s(60)=c(11)*d(12)+c(12)*d(11)+c(14)*d(16)+c(16)*d(14)                     
      s(61)=c(11)*d(13)+c(13)*d(11)+c(15)*d(18)+c(18)*d(15)                     
      s(62)=c(12)*d(13)+c(13)*d(12)+c(17)*d(19)+c(19)*d(17)                     
      s(63)=c(11)*d(14)+c(14)*d(11)                                             
      s(64)=c(11)*d(15)+c(15)*d(11)                                             
      s(65)=c(13)*d(18)+c(18)*d(13)                                             
      s(66)=c(13)*d(19)+c(19)*d(13)                                             
      s(67)=c(12)*d(17)+d(12)*c(17)                                             
      s(68)=c(12)*d(16)+c(16)*d(12)                                             
      s(69)=c(11)*d(16)+c(16)*d(11)+c(14)*d(14)                                 
      s(70)=c(11)*d(18)+c(18)*d(11)+c(15)*d(15)                                 
      s(71)=c(13)*d(15)+c(15)*d(13)+c(18)*d(18)                                 
      s(72)=c(13)*d(17)+c(17)*d(13)+c(19)*d(19)                                 
      s(73)=c(12)*d(14)+c(14)*d(12)+c(16)*d(16)                                 
      s(74)=c(12)*d(19)+c(19)*d(12)+c(17)*d(17)                                 
      s(75)=c(11)*d(17)+c(17)*d(11)+c(14)*d(20)+c(20)*d(14)+                    
     1      c(15)*d(16)+c(16)*d(15)                                             
      s(76)=c(11)*d(19)+c(19)*d(11)+c(20)*d(15)+c(15)*d(20)+                    
     1      c(14)*d(18)+c(18)*d(14)                                             
      s(77)=c(12)*d(18)+c(18)*d(12)+c(16)*d(19)+c(19)*d(16)+                    
     1      c(17)*d(20)+c(20)*d(17)                                             
      s(78)=c(13)*d(14)+c(14)*d(13)+c(15)*d(19)+c(19)*d(15)+                    
     1      c(18)*d(20)+c(20)*d(18)                                             
      s(79)=c(12)*d(15)+c(15)*d(12)+c(14)*d(17)+c(17)*d(14)+                    
     1      c(16)*d(20)+c(20)*d(16)                                             
      s(80)=c(13)*d(16)+c(16)*d(13)+c(17)*d(18)+c(18)*d(17)+                    
     1      c(19)*d(20)+c(20)*d(19)                                             
      s(81)=c(11)*d(20)+c(20)*d(11)+c(14)*d(15)+c(15)*d(14)                     
      s(82)=c(12)*d(20)+c(20)*d(12)+c(16)*d(17)+c(17)*d(16)                     
      s(83)=c(13)*d(20)+c(20)*d(13)+c(18)*d(19)+c(19)*d(18)                     
      s(84)=c(14)*d(19)+c(19)*d(14)+c(15)*d(17)+c(17)*d(15)+                    
     1      c(16)*d(18)+c(18)*d(16)+c(20)*d(20)                                 
c              end of f - f                                                     
      return                                                                    
      end                                                                       

      subroutine lmnpre(l,m,n,lmnexp,lmnfak)                          
      implicit real*8(a-h,o-z)                                                  
      real*8  lmnfak
      integer lmnexp
      dimension dftr(7)                                                         
      data dftr /1.d0,1.d0,3.d0,15.d0,105.d0,945.d0,10395.d0/                   
      lmnfak=0
      lmnexp=0
      lh=l/2                                                                    
      if(2*lh-l) 50,1,50                                                        
   1  mh=m/2                                                                    
      if(2*mh-m) 50,2,50                                                        
   2  nh=n/2                                                                    
      if(2*nh-n) 50,3,50                                                        
   3  lmnexp=lh+mh+nh                      
      lmnfak=dftr(lh+1)*dftr(mh+1)*dftr(nh+1)                                    
  50  return                                                                    
      end                                                                       

      double precision function olap2(l,m,n,arg,gama)                                a
      implicit real*8(a-h,o-z)                                                  
      dimension dftr(7)                                                         
      data dftr /1.d0,1.d0,3.d0,15.d0,105.d0,945.d0,10395.d0/                   
      lh=l/2                                                                    
      if(2*lh-l) 50,1,50                                                        
   1  mh=m/2                                                                    
      if(2*mh-m) 50,2,50                                                        
   2  nh=n/2                                                                    
      if(2*nh-n) 50,3,50                                                        
   3  olap2=arg*gama**(lh+mh+nh)*                       
     1      dftr(lh+1)*dftr(mh+1)*dftr(nh+1)                                    
      return                                                                    
  50  olap2=0.d0                                                                 
      return                                                                    
      end                                                                       

      double precision function olap(l,m,n,gama)                                
      implicit real*8(a-h,o-z)                                                  
      dimension dftr(7)                                                         
      data pi/3.1415926535897932384626433832795029d0/
      data dftr /1.d0,1.d0,3.d0,15.d0,105.d0,945.d0,10395.d0/                   
      data thehlf /0.5d0/                                                       
      gm=1.0d0/gama
      lh=l/2                                                                    
      if(2*lh-l) 50,1,50                                                        
   1  mh=m/2                                                                    
      if(2*mh-m) 50,2,50                                                        
   2  nh=n/2                                                                    
      if(2*nh-n) 50,3,50                                                        
   3  olap=(dsqrt(pi*gm))**3*(thehlf*gm)**(lh+mh+nh)*                       
     1      dftr(lh+1)*dftr(mh+1)*dftr(nh+1)                                    
      return                                                                    
  50  olap=0.d0                                                                 
      return                                                                    
      end                                                                       

      subroutine opab1(l,m,n,ga,v,d)
c           electronic part of dipole moment
      implicit real*8(a-h,o-z)
      dimension v(*),d(*),g(4)
      g(1)=olap(l,m,n,ga)
      g(2)=olap(l+1,m,n,ga)
      g(3)=olap(l,m+1,n,ga)
      g(4)=olap(l,m,n+1,ga)
      v(1)=g(2)+d(1)*g(1)
      v(2)=g(3)+d(2)*g(1)
      v(3)=g(4)+d(3)*g(1)
c     do 300 i=1,3         CSG sign moved to calling routine
c 300 v(i)=-v(i)
      return
      end

      subroutine opab4(l,m,n,ga,v,d)
c           electronic part of second moment, x^2,y^2,z^2 only
      implicit real*8(a-h,o-z)
      dimension v(*),d(*),g(10)
      g(1)=olap(l,m,n,ga)
      g(2)=olap(l+1,m,n,ga)
      g(3)=olap(l,m+1,n,ga)
      g(4)=olap(l,m,n+1,ga)
      g(5)=olap(l+2,m,n,ga)
      g(6)=olap(l,m+2,n,ga)
      g(7)=olap(l,m,n+2,ga)
!     g(8)=olap(l+1,m+1,n,ga)
!     g(9)=olap(l+1,m,n+1,ga)
!     g(10)=olap(l,m+1,n+1,ga)
      v(1)=g(5)+2.d0*d(1)*g(2)+d(1)**2*g(1)
      v(2)=g(6)+2.d0*d(2)*g(3)+d(2)**2*g(1)
      v(3)=g(7)+2.d0*d(3)*g(4)+d(3)**2*g(1)
!     v(4)=g(8)+d(1)*g(3)+d(2)*g(2)+d(1)*d(2)*g(1)
!     v(5)=g(9)+d(1)*g(4)+d(3)*g(2)+d(1)*d(3)*g(1)
!     v(6)=g(10)+d(2)*g(4)+d(3)*g(3)+d(2)*d(3)*g(1)
!     do 300 i=1,3
!300  v(i)=-v(i)
      return
      end

      subroutine opab46(l,m,n,ga,v,d)
c           electronic part of second moment, x^2,y^2,z^2 only
      implicit real*8(a-h,o-z)
      dimension v(*),d(*),g(10)
      g(1)=olap(l,m,n,ga)
      g(2)=olap(l+1,m,n,ga)
      g(3)=olap(l,m+1,n,ga)
      g(4)=olap(l,m,n+1,ga)
      g(5)=olap(l+2,m,n,ga)
      g(6)=olap(l,m+2,n,ga)
      g(7)=olap(l,m,n+2,ga)
      g(8)=olap(l+1,m+1,n,ga)
      g(9)=olap(l+1,m,n+1,ga)
      g(10)=olap(l,m+1,n+1,ga)
      v(4)=g(5)+2.d0*d(1)*g(2)+d(1)**2*g(1)
      v(5)=g(6)+2.d0*d(2)*g(3)+d(2)**2*g(1)
      v(6)=g(7)+2.d0*d(3)*g(4)+d(3)**2*g(1)
      v(7)=g(8)+d(1)*g(3)+d(2)*g(2)+d(1)*d(2)*g(1)
      v(8)=g(9)+d(1)*g(4)+d(3)*g(2)+d(1)*d(3)*g(1)
      v(9)=g(10)+d(2)*g(4)+d(3)*g(3)+d(2)*d(3)*g(1)
      return
      end

      subroutine opac3(l,m,n,ga,v,d)
c           electronic part of charge density
      implicit real*8(a-h,o-z)
      dimension v(1),d(1),dd(3)
      ddsq=0.d0
      do 1 i=1,3
      dd(i)=-d(i)
      ddsq=ddsq+d(i)**2
    1 continue
      v(1) =1.d0
      if(l) 2,3,2
    2 v(1)=dd(1)**l
    3 if(m) 4,5,4
    4 v(1)=v(1)*dd(2)**m
    5 if(n) 6,7,6
    6 v(1)=v(1)*dd(3)**n
    7 v(1)=v(1)*dexp(-ga*ddsq)
      return
      end

      subroutine opad1(l,m,n,ga,v,d)                                            
c           electronic part of overlap                                          
      implicit real*8(a-h,o-z)                                                  
      dimension v(1),d(1)                                                       
      v(1)=olap(l,m,n,ga)                                                       
      return                                                                    
      end                                                                       
                                                                                
      subroutine opap4(l,m,n,gc,v,rc)
c
c     this subroutine calculates the expectation-values of
c     -p**2/2 (kinetic energy) and p**4/(8*cl**2) (relativistic
c     correction of the kinetic energy)
c
c     the program works correct only for s and p and the irreducible
c     linear combinations of the primitive d,f,... functions
c
c     written by s. brode in april,1984
c
      implicit real*8 (a-h,o-z)
      logical lodd,modd,nodd,leven,meven,neven
      dimension v(1),rc(3),rca(3),rcb(3),srcab(3),
     &lmni(35),o(25),dftr(7)
c ----- ra,rb centres ga,gb exponents ia,ib monom indices
      common /abfunc/ ra(3),rb(3),ga,gb,ia,ib
      data a0/0.0d0/,a1/1.0d0/,a2/2.0d0/,a4/4.0d0/,a8/8.0d0/,
     &lmni/0,3*1,6*2,10*3,15*4/,cl/137.03604d0/,
     &dftr/1.0d0,1.0d0,3.0d0,15.0d0,105.0d0,945.0d0,103935.0d0/
      data pi/3.1415926535897932384626433832795029d0/
c
c     define the overlap-integral (in line function definition)
c
      ola(ll,mm,nn)=gch**(ll/2+mm/2+nn/2)*dftr(ll/2+1)*dftr(mm/2+1)*
     &dftr(nn/2+1)
c
c     some previous calculations
c
c     do 10 i=1,25
c  10 o(i)=a0
      o(1:7)=0
      gch=a1/(a2*gc)
      ropigc=dsqrt(pi/gc)*pi/gc
      lodd=mod(l,2).eq.1
      leven=.not.lodd
      modd=mod(m,2).eq.1
      meven=.not.modd
      nodd=mod(n,2).eq.1
      neven=.not.nodd
      ga2=ga*a2
      gb2=gb*a2
      gab4=ga2*gb2
      prcaa=a0
      prcbb=a0
      do 20 i=1,3
      rca(i)=rc(i)-ra(i)
      rcb(i)=rc(i)-rb(i)
      srcab(i)=rca(i)+rcb(i)
      prcaa=prcaa+rca(i)*rca(i)
      prcbb=prcbb+rcb(i)*rcb(i)
   20 continue
      zwa=ga2*prcaa-dble(2*lmni(ia)+3)
      zwb=gb2*prcbb-dble(2*lmni(ib)+3)
c
c     calculation of the overlap-integrals
c
      if(lodd.or.modd.or.nodd) go to 110
      o( 1)=ola(l   ,m  ,n  )
      o( 5)=ola(l+2,m   ,n  )
      o( 6)=ola(l   ,m+2,n  )
      o( 7)=ola(l   ,m  ,n+2)
c     o(11)=ola(l+4,m   ,n  )
c     o(12)=ola(l   ,m+4,n  )
c     o(13)=ola(l   ,m  ,n+4)
c     o(20)=ola(l+2,m+2,n   )
c     o(21)=ola(l+2,m   ,n+2)
c     o(22)=ola(l   ,m+2,n+2)
  110 continue
c
      if(leven.or.modd.or.nodd) go to 120
      o( 2)=ola(l+1,m   ,n  )
c     o( 8)=ola(l+3,m   ,n  )
c     o(14)=ola(l+1,m+2,n   )
c     o(15)=ola(l+1,m   ,n+2)
  120 continue
c
      if(lodd.or.meven.or.nodd) go to 130
      o( 3)=ola(l   ,m+1,n  )
c     o( 9)=ola(l   ,m+3,n  )
c     o(16)=ola(l+2,m+1,n   )
c     o(17)=ola(l   ,m+1,n+2)
  130 continue
c
      if(lodd.or.modd.or.neven) go to 140
      o( 4)=ola(l   ,m  ,n+1)
c     o(10)=ola(l   ,m  ,n+3)
c     o(18)=ola(l+2,m   ,n+1)
c     o(19)=ola(l   ,m+2,n+1)
  140 continue
c
c     if(leven.or.meven.or.nodd) go to 150
c     o(23)=ola(l+1,m+1,n   )
c 150 continue
c     if(leven.or.modd.or.neven) go to 160
c     o(24)=ola(l+1,m   ,n+1)
c 160 continue
c     if(lodd.or.meven.or.neven) go to 170
c     o(25)=ola(l   ,m+1,n+1)
c 170 continue
c
c     calculation of kinetic energy
c
      p2=-ga*(zwa*o(1)
     &       +ga2*(a2*(rca(1)*o(2)+rca(2)*o(3)+rca(3)*o(4))
     &       +o(5)+o(6)+o(7)))

      
c     calculation of relativistic correction
 
c     p4=(gab4*(o(11)+o(12)+o(13)+a2*(o(20)+o(21)+o(22))
c    &         +a2*(srcab(1)*(o( 8)+o(16)+o(18))
c    &             +srcab(2)*(o(14)+o( 9)+o(19))
c    &             +srcab(3)*(o(15)+o(17)+o(10)))
c    &         +a4*(rca(1)*rcb(1)*o(5)
c    &             +rca(2)*rcb(2)*o(6)
c    &             +rca(3)*rcb(3)*o(7)
c    &             +a2*((rca(1)*rcb(2)+rca(2)+rcb(1))*o(23)
c    &                 +(rca(1)*rcb(3)+rca(3)+rcb(1))*o(24)
c    &                 +(rca(2)*rcb(3)+rca(3)+rcb(2))*o(25))))
c    &   +(gb2*zwa+ga2*zwb)*(o(5)+o(6)+o(7))
c    &   +a2*(gb2*zwa*(rcb(1)*o(2)+rcb(2)*o(3)+rcb(3)*o(4))
c    &       +ga2*zwb*(rca(1)*o(2)+rca(2)*o(3)+rca(3)*o(4)))
c    &   +zwa*zwb*o(1))
c    &   *(-gab4)/(a8*cl**2)
 
c     store p2 and p4 to the correct storage-locations
 
      v(1)=p2*ropigc
c     v(2)=p4*ropigc
c     v(3)=v(1)+v(2)
      return
      end

      subroutine opaa0(l,m,n,ga,v,d)
c           electronic part of potential
      implicit real*8(a-h,o-z)
      dimension v(*),d(*),fnu(7),fn(7),fd(7),dp(3)
      data fn /1.d0,1.d0,2.d0,6.d0,24.d0,120.d0,
     1 720.d0/
      data fd/1.d0,1.d0,0.5d0,0.1666666666667d0,0.416666666667d-1,
     1 .833333333333d-2,.1388888888889d-2/
      data pi/3.1415926535897932384626433832795029d0/
      itot=l+m+n
      itoth=itot/2
      pre=(2.d0*pi/ga)*fn(l+1)*fn(m+1)*fn(n+1)
      if(2*itoth-itot) 1,2,1
   1  pre=-pre
   2  del=.25d0/ga
      x=ga*(d(1)**2+d(2)**2+d(3)**2)
      xx=2.d0*x
      call fmc(6,x,expmx,fmch)
      fnu(7)=fmch
      fnu(6)=(expmx+xx*fnu(7))/11.d0
      fnu(5)=(expmx+xx*fnu(6))/9.d0
      fnu(4)=(expmx+xx*fnu(5))/7.d0
      fnu(3)=(expmx+xx*fnu(4))/5.d0
      fnu(2)=(expmx+xx*fnu(3))/3.d0
      fnu(1)=expmx+xx*fnu(2)
      dp(1)=1.d0
      dp(2)=del
      dp(3)=del**2
      v(1)= - pre*aainer(0,0,0,l,m,n,d,dp,fnu,fn,fd)
      return
      end

      double precision function aainer(r,s,t,l,m,n,d,dp,fnu,fn,fd)
      implicit real*8(a-h,o-z)
      dimension  d(*),dp(*),fnu(*),fn(*),fd(*)
      integer r,s,t,u,v,w,uvwt,rstt,uvwth
      rstt=r+s+t
      lmnt=l+m+n
      lmnrst=rstt+lmnt
      lh=l/2
      mh=m/2
      nh=n/2
      aainer=0.0d0
      last1 = lh+1
      do 6 ii=1,last1
       i = ii-1
       il=l-2*i+1
       ilr=r+il
       ilrh=(ilr-1)/2
!     if (ilr.gt.9 .or. i+1.gt.9 .or. il.gt.9) stop 'aainer: ilr'
       fi=fn(ilr)*fd(il)*fd(i+1)
       last2 = mh+1
       do 5 jj=1,last2
        j = jj-1
        jm=m-2*j+1
        jms=s+jm
        jmsh=(jms-1)/2
!     if (jms.gt.9 .or. j+1.gt.9 .or. jm.gt.9) stop 'aainer: ilr'
        fij=fn(jms)*fd(jm)*fd(j+1)*fi
        last3 = nh+1
        do 4 kk=1,last3
         k = kk-1
         kn=n-2*k+1
         knt=t+kn
         knth=(knt-1)/2
         ijkt=i+j+k
!     if (knt.gt.9 .or. k+1.gt.9 .or. kn.gt.9) stop 'aainer: ilr'
         fijk=fn(knt)*fd(kn)*fd(k+1)*dp(ijkt+1)*fij
         lmrsij=lmnrst-2*ijkt
         last4 = ilrh+1
         do 3 iu=1,last4
          u = iu-1
          ilru=ilr-2*u
!     if (u+1.gt.9 .or. ilru.gt.9) stop 'aainer: ilru'
          fu=fd(u+1)*fd(ilru)
          if(dabs(d(1))-0.000001d0)10,10,11
   10     if(ilru-1) 12,12, 3
   11     fu=fu*d(1)**(ilru-1)
   12     last5 = jmsh+1
          do 2 iv=1,last5
           v = iv-1
           jmsv=jms-2*v
           fuv=fu*fd(v+1)*fd(jmsv)
           if(dabs(d(2))-0.000001d0)20,20,21
   20      if(jmsv-1) 22,22,2
   21      fuv=fuv*d(2)**(jmsv-1)
   22      last6 = knth+1
           do 1 iw=1,last6
            w = iw-1
            kntw=knt-2*w
            fuvw=fuv*fd(w+1) *fd(kntw)
            if(dabs(d(3))-0.000001d0)30,30,31
   30       if(kntw-1) 32,32,1
   31       fuvw=fuvw*d(3)**(kntw-1)
   32       uvwt=u+v+w
            uvwth=uvwt/2
            if(2*uvwth-uvwt) 33,34,33
   33       fuvw=-fuvw
   34       nuindx=lmrsij-uvwt
            fuvw=fuvw*fnu(nuindx  +1)*dp(uvwt+1)
            aainer=fijk*fuvw+aainer
    1      continue
    2     continue
    3    continue
    4   continue
    5  continue
    6 continue
      return
      end

      subroutine fmc(mvar,xvar,expmx,fmch)
      implicit real*8(a-h,o-z)
      m=mvar
      x=xvar
cmk ---- prevent underflow error messages
      expmx=0.d0
      if(x.lt.50.d0) expmx=dexp(-x)
c     expmx=dexp(-x)
c
      if(x -20.d0) 10,10,20
c         x less than or equal to 20.0.
   10 a=m
      a=a+0.5d0
      term=1.d0/a
      ptlsum=term
      do 11 i=2,60
       a=a+1.d0
       term=term*x/a
       ptlsum=ptlsum+term
       if (term/ptlsum- 1.d-8    ) 12, 11, 11
   11  continue
       print 999,m,x
  12  fmch=.5d0*ptlsum*expmx
      return
c         x greater than 20.0.
   20 a=m
      xd=1.d0/x
      b=a+0.5d0
      a=a-0.5d0
      approx=.886226925428d0*(dsqrt(xd)*xd**m)
      if (m) 21, 23, 21
   21 do 22 i=1,m
       b=b-1.d0
   22 approx=approx*b
  23  fimult=.5d0*expmx*xd
      fiprop=fimult/approx
      term=1.d0
      ptlsum=term
      notrms=x
      notrms=notrms+m
      do 24 i=2,notrms
       term=term*a*xd
       ptlsum=ptlsum+term
      if (dabs(term*fiprop/ptlsum)-1.d-8    ) 25,25,24
   24 a=a-1.d0
      print 999,m,x
   25 fmch=approx-fimult*ptlsum
      return
  999 format (24h no convergence for fmch, i6, d16.9)
      end

