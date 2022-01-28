subroutine cma(n,at,coord,pnt)        
      use com      ! general stuff         
      implicit none
      integer n,at(n)
      real*8  coord(3,n),pnt(3)

      real*8 sumw,sumwx,sumwy,sumwz,atmass
      integer i
                   
      sumw =1.d-20   
      sumwx=0.d0  
      sumwy=0.d0 
      sumwz=0.d0

      do i=1,n                                                        
         atmass=ams(at(i))                                                  
         sumw=sumw+atmass                                                    
         sumwx=sumwx+atmass*coord(1,i)                                       
         sumwy=sumwy+atmass*coord(2,i)
         sumwz=sumwz+atmass*coord(3,i)
      enddo                                                                     
                                                                                
      pnt(1)=sumwx/sumw                                                          
      pnt(2)=sumwy/sumw                                                          
      pnt(3)=sumwz/sumw                                                          
end                                                                            

