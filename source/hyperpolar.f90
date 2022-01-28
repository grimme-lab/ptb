      call calcrab(n,at,xyz,rab)
      call modbas(n,at,4) 
      call sint(n,ndim,at,xyz,rab,S,xnorm)       ! exact S
      y=0.0002_wp ! field strength
      do m=1,3
         efield=0_wp           
         efield(m)=y           
         call shellq(.false.,2,n,ndim,nel,nopen,ihomo,at,chrg,xyz,z,rab,pnt,xnorm,S,&
     &               efield,psh,q,P,wbo,dip,alpr)
         efield(m)=-y          
         call shellq(.false.,2,n,ndim,nel,nopen,ihomo,at,chrg,xyz,z,rab,pnt,xnorm,S,&
     &               efield,psh,q,P,wbo,dip,alp)
         efield=0_wp           
         beta(1:6,m)=(alpr(1:6)-alp(1:6))/(2_wp*y)
      enddo
      
      k=0
      do i=1,n
         do j=1,3
            k=k+1
            y=0.0001_wp  ! field strength
            fdgrad(1:3,k,7)=(fdgrad(1:3,k,1)-fdgrad(1:3,k,2))/(2_wp*y)
            fdgrad(1:3,k,8)=(fdgrad(1:3,k,3)-fdgrad(1:3,k,4))/(2_wp*y)
            fdgrad(1:3,k,9)=(fdgrad(1:3,k,5)-fdgrad(1:3,k,6))/(2_wp*y)
            y=0.0020_wp ! field strength TM script
            fdgrad_ref(1:3,k,7)=(fdgrad_ref(1:3,k,1)-fdgrad_ref(1:3,k,2))/(2_wp*y)
            fdgrad_ref(1:3,k,8)=(fdgrad_ref(1:3,k,3)-fdgrad_ref(1:3,k,4))/(2_wp*y)
            fdgrad_ref(1:3,k,9)=(fdgrad_ref(1:3,k,5)-fdgrad_ref(1:3,k,6))/(2_wp*y)
         enddo
      enddo
