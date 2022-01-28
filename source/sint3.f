      SUBROUTINE ss1e(xyza,xyzb,p,ab,s00,D)
      IMPLICIT REAL*8(A-H,O-Z)                            
      dimension d(6,6),xyza(3),xyzb(3),p(3)

      d = 0
      d(1,1)=s00

      end   

      SUBROUTINE sp1e(xyza,xyzb,p,ab,s00,D)
      IMPLICIT REAL*8(A-H,O-Z)                            
      dimension d(6,6),xyza(3),xyzb(3),p(3)

      d = 0
      d(1,1)=(p(1)-xyzb(1))*s00
      d(2,1)=(p(2)-xyzb(2))*s00
      d(3,1)=(p(3)-xyzb(3))*s00

      end   

      SUBROUTINE ps1e(xyza,xyzb,p,ab,s00,D)
      IMPLICIT REAL*8(A-H,O-Z)                            
      dimension d(6,6),xyza(3),xyzb(3),p(3)

      d = 0
      d(1,1)=(p(1)-xyza(1))*s00
      d(1,2)=(p(2)-xyza(2))*s00
      d(1,3)=(p(3)-xyza(3))*s00

      end   

      SUBROUTINE sd1e(xyza,xyzb,p,ab,s00,D)
      IMPLICIT REAL*8(A-H,O-Z)                            
      dimension d(6,6),xyza(3),xyzb(3),p(3)

      d = 0
      dum=0.5*ab*s00
      dp1=p(1)-xyzb(1)
      dp2=p(2)-xyzb(2)
      dp3=p(3)-xyzb(3)
      sp1 =s00*dp1            
      sp2 =s00*dp2           
      sp3 =s00*dp3           
      d(1,1)=    dp1       *sp1+dum       
      d(2,1)=    dp2       *sp2+dum        
      d(3,1)=    dp3       *sp3+dum       
      d(4,1)=    dp2       *sp1
      d(5,1)=    dp3       *sp1
      d(6,1)=    dp3       *sp2

      end   

      SUBROUTINE ds1e(xyza,xyzb,p,ab,s00,D)
      IMPLICIT REAL*8(A-H,O-Z)                            
      dimension d(6,6),xyza(3),xyzb(3),p(3)

      d = 0
      dum=0.5*ab*s00
      dp1=p(1)-xyza(1)
      dp2=p(2)-xyza(2)
      dp3=p(3)-xyza(3)
      sp1 =s00*dp1            
      sp2 =s00*dp2           
      sp3 =s00*dp3           
      d(1,1)=    dp1       *sp1+dum       
      d(1,2)=    dp2       *sp2+dum        
      d(1,3)=    dp3       *sp3+dum       
      d(1,4)=    dp2       *sp1
      d(1,5)=    dp3       *sp1
      d(1,6)=    dp3       *sp2

      end   

      SUBROUTINE pp1e(xyza,xyzb,p,ab,s00,D)
      IMPLICIT REAL*8(A-H,O-Z)                            
      dimension d(6,6),xyza(3),xyzb(3),p(3)

      d = 0
      ab05=ab*0.5
      do id = 1, 3
         xps=(p(id)-xyza(id))*s00         
         d(1,id)=(p(1)-xyzb(1))*xps
         d(2,id)=(p(2)-xyzb(2))*xps
         d(3,id)=(p(3)-xyzb(3))*xps
         d(id,id)=d(id,id)+ab05*s00
      enddo

      end   

      SUBROUTINE pd1e(xyza,xyzb,p,ab,s00,D)
      IMPLICIT REAL*8(A-H,O-Z)                            
      dimension d(6,6),xyza(3),xyzb(3),p(3)

      d = 0
      ab05=ab*0.5
      dum=ab05*s00
      dp1=p(1)-xyzb(1)
      dp2=p(2)-xyzb(2)
      dp3=p(3)-xyzb(3)
      sp1=s00*dp1
      sp2=s00*dp2
      sp3=s00*dp3
      sdxx=    dp1       *sp1+dum       
      sdyy=    dp2       *sp2+dum        
      sdzz=    dp3       *sp3+dum       
      sdxy=    dp2       *sp1
      sdxz=    dp3       *sp1
      sdyz=    dp3       *sp2

      dp1=p(1)-xyza(1)
      d(1,1)=dp1*sdxx+  ab*sp1 
      d(2,1)=dp1*sdyy
      d(3,1)=dp1*sdzz
      d(4,1)=dp1*sdxy+ab05*sp2
      d(5,1)=dp1*sdxz+ab05*sp3
      d(6,1)=dp1*sdyz

      dp1=p(2)-xyza(2)
      d(1,2)=dp1*sdxx
      d(2,2)=dp1*sdyy+  ab*sp2
      d(3,2)=dp1*sdzz
      d(4,2)=dp1*sdxy+ab05*sp1
      d(5,2)=dp1*sdxz
      d(6,2)=dp1*sdyz+ab05*sp3

      dp1=p(3)-xyza(3)
      d(1,3)=dp1*sdxx
      d(2,3)=dp1*sdyy
      d(3,3)=dp1*sdzz+  ab*sp3
      d(4,3)=dp1*sdxy
      d(5,3)=dp1*sdxz+ab05*sp1
      d(6,3)=dp1*sdyz+ab05*sp2

      end   

      SUBROUTINE dp1e(xyza,xyzb,p,ab,s00,D)
      IMPLICIT REAL*8(A-H,O-Z)                            
      dimension d(6,6),xyza(3),xyzb(3),p(3)

      d = 0
      ab05=ab*0.5
      dum=ab05*s00
      dp1=p(1)-xyza(1)
      dp2=p(2)-xyza(2)
      dp3=p(3)-xyza(3)
      sp1=s00*dp1
      sp2=s00*dp2
      sp3=s00*dp3
      sdxx=    dp1       *sp1+dum       
      sdyy=    dp2       *sp2+dum        
      sdzz=    dp3       *sp3+dum       
      sdxy=    dp2       *sp1
      sdxz=    dp3       *sp1
      sdyz=    dp3       *sp2

      dp1=p(1)-xyzb(1)
      d(1,1)=dp1*sdxx+  ab*sp1 
      d(1,2)=dp1*sdyy
      d(1,3)=dp1*sdzz
      d(1,4)=dp1*sdxy+ab05*sp2
      d(1,5)=dp1*sdxz+ab05*sp3
      d(1,6)=dp1*sdyz

      dp1=p(2)-xyzb(2)
      d(2,1)=dp1*sdxx
      d(2,2)=dp1*sdyy+  ab*sp2
      d(2,3)=dp1*sdzz
      d(2,4)=dp1*sdxy+ab05*sp1
      d(2,5)=dp1*sdxz
      d(2,6)=dp1*sdyz+ab05*sp3

      dp1=p(3)-xyzb(3)
      d(3,1)=dp1*sdxx
      d(3,2)=dp1*sdyy
      d(3,3)=dp1*sdzz+  ab*sp3
      d(3,4)=dp1*sdxy
      d(3,5)=dp1*sdxz+ab05*sp1
      d(3,6)=dp1*sdyz+ab05*sp2

      end   

      SUBROUTINE dd1e(xyza,xyzb,p,ab,s00,D)
      IMPLICIT REAL*8(A-H,O-Z)                            
      dimension d(6,6),xyza(3),xyzb(3),p(3)
      dimension pip(3),pjp(3),deltapb(3),ps(3),dp(3)
      dimension ii(6),jj(6)
      data ii/1,2,3,1,1,2/
      data jj/1,2,3,2,3,3/

      d = 0
      ab05=ab*0.5

      do ipi=1,6
         iii=ii(ipi)
         jjj=jj(ipi)

         fij=0.0d0
         if(iii.eq.jjj) fij=ab05
         ds=(p(jjj)-xyza(jjj))*(p(iii)-xyza(iii))*s00+fij*s00
         ps(1)=(p(1)-xyza(1))*s00
         ps(2)=(p(2)-xyza(2))*s00
         ps(3)=(p(3)-xyza(3))*s00
         deltapb(1)=(p(1)-xyzb(1))
         deltapb(2)=(p(2)-xyzb(2))
         deltapb(3)=(p(3)-xyzb(3))
         do m=1,3
            fik=0.0d0
            fjk=0.0d0
            if(iii.eq.m) fik=ab05
            if(jjj.eq.m) fjk=ab05  
            dp(m)=deltapb(m)*ds+fik*ps(jjj)+fjk*ps(iii)  
         enddo
         do m=1,3
            pip(m)=deltapb(m)*ps(iii)
            if(iii.eq.m) pip(m)=pip(m)+ab05*s00
         enddo
         if(iii.ne.jjj) then
            do m=1,3
               pjp(m)=deltapb(m)*ps(jjj)
               if(jjj.eq.m) pjp(m)=pjp(m)+ab05*s00
            enddo
         else
            pjp(1)=pip(1)
            pjp(2)=pip(2)
            pjp(3)=pip(3)
         endif
         ab05ds=ab05*ds
         d(1,ipi)=deltapb(1)    *dp(1)+ab05ds 
         d(2,ipi)=deltapb(2)    *dp(2)+ab05ds 
         d(3,ipi)=deltapb(3)    *dp(3)+ab05ds
         d(4,ipi)=deltapb(2)    *dp(1) 
         d(5,ipi)=deltapb(3)    *dp(1) 
         d(6,ipi)=deltapb(3)    *dp(2) 
         do m=1,6
            kkk=ii(m)
            lll=jj(m)
            if(iii.eq.lll) d(m,ipi)=d(m,ipi)+ab05*pjp(kkk)
            if(jjj.eq.lll) d(m,ipi)=d(m,ipi)+ab05*pip(kkk)
         enddo

      enddo

      end
